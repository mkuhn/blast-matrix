/* $Id: tar.cpp 257864 2011-03-16 15:10:08Z rafanovi $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vladimir Ivanov
 *           Anton Lavrentiev
 *
 * File Description:
 *   Tar archive API.
 *
 *   Supports subsets of POSIX.1-1988 (ustar), POSIX 1003.1-2001 (posix),
 *   old GNU (POSIX 1003.1), and V7 formats (all partially but reasonably).
 *   New archives are created using POSIX (genuine ustar) format, using
 *   GNU extensions for long names/links only when unavoidable.
 *   Cannot handle the exotics like sparse / contiguous files,
 *   multivolume / incremental archives, etc, but just regular files,
 *   devices (character or block), FIFOs, directories, and limited links:
 *   can extract both hard- and symlinks, but can store symlinks only.
 *   This version is minimally PAX (Partable Archive Interchange) aware
 *   for file extractions (but cannot use PAX extensions to store files).
 *
 */

#include <ncbi_pch.hpp>
/* Cancel __wur (warn unused result) ill effects in GCC */
#ifdef   _FORTIFY_SOURCE
#  undef _FORTIFY_SOURCE
#endif /*_FORTIFY_SOURCE*/
#define  _FORTIFY_SOURCE 0
#include <corelib/ncbi_limits.h>
#include <corelib/ncbimisc.hpp>
#include <corelib/ncbi_system.hpp>
#include <util/compress/tar.hpp>
#include <util/error_codes.hpp>

#if !defined(NCBI_OS_MSWIN)  &&  !defined(NCBI_OS_UNIX)
#  error "Class CTar can be defined on MS-Windows and UNIX platforms only!"
#endif

#if   defined(NCBI_OS_UNIX)
#  include <unistd.h>
#  ifdef NCBI_OS_IRIX
#    include <sys/mkdev.h>
#    if !defined(major)  ||  !defined(minor)  ||  !defined(makedev)
#      error "Device macros undefined in this UNIX build!"
#    endif
#  endif
#elif defined(NCBI_OS_MSWIN)
#  include <io.h>
typedef unsigned int mode_t;
typedef short        uid_t;
typedef short        gid_t;
#endif // NCBI_OS


#define NCBI_USE_ERRCODE_X  Util_Compress
#define NCBI_MODULE         NCBITAR


BEGIN_NCBI_SCOPE


/////////////////////////////////////////////////////////////////////////////
//
// TAR helper routines
//

// Convert a number to an octal string padded to the left
// with [leading] zeros ('0') and having _no_ trailing '\0'.
static bool s_NumToOctal(Uint8 val, char* ptr, size_t len)
{
    _ASSERT(len > 0);
    do {
        ptr[--len] = '0' + char(val & 7);
        val >>= 3;
    } while (len);
    return val ? false : true;
}


// Convert an octal number (possibly preceded by spaces) to numeric form.
// Stop either at the end of the field or at first '\0' (if any).
static bool s_OctalToNum(Uint8& val, const char* ptr, size_t len)
{
    _ASSERT(ptr  &&  len > 0);
    size_t i = *ptr ? 0 : 1;
    while (i < len  &&  ptr[i]) {
        if (!isspace((unsigned char) ptr[i]))
            break;
        i++;
    }
    val = 0;
    bool okay = false;
    while (i < len  &&  ptr[i] >= '0'  &&  ptr[i] <= '7') {
        okay  = true;
        val <<= 3;
        val  |= ptr[i++] - '0';
    }
    while (i < len  &&  ptr[i]) {
        if (!isspace((unsigned char) ptr[i]))
            return false;
        i++;
    }
    return okay;
}


static bool s_NumToBase256(Uint8 val, char* ptr, size_t len)
{
    _ASSERT(len > 0);
    do {
        ptr[--len] = (unsigned char)(val & 0xFF);
        val >>= 8;
    } while (len);
    *ptr |= '\x80';  // set base-256 encoding flag
    return val ? false : true;
}


// Return 0 (false) if conversion failed; 1 if the value converted to
// conventional octal representation (perhaps, with terminating '\0'
// sacrificed), or -1 if the value converted using base-256.
static int s_EncodeUint8(Uint8 val, char* ptr, size_t len)
{                                          // Max file size:
    if (s_NumToOctal  (val, ptr,   len)) { //   8GiB-1
        return  1/*okay*/;
    }
    if (s_NumToOctal  (val, ptr, ++len)) { //   64GiB-1
        return  1/*okay*/;
    }
    if (s_NumToBase256(val, ptr,   len)) { //   up to 2^94-1
        return -1/*okay, base-256*/;
    }
    return 0/*failure*/;
}


// Return true if conversion succeeded;  false otherwise.
static bool s_Base256ToNum(Uint8& val, const char* ptr, size_t len)
{
    const Uint8 lim = kMax_UI8 >> 8;
    if (*ptr & '\x40') { /*negative base-256?*/
        return false;
    }
    val = *ptr++ & '\x3F';
    while (--len) {
        if (val > lim) {
            return false;
        }
        val <<= 8;
        val  |= (unsigned char)(*ptr++);
    }
    return true;
}


// Return 0 (false) if conversion failed; 1 if the value was read into
// as a conventional octal string (perhaps, without the terminating '\0');
// or -1 if base-256 representation used.
static int s_DecodeUint8(Uint8& val, const char* ptr, size_t len)
{
    if (*ptr & '\x80') {
        return s_Base256ToNum(val, ptr, len) ? -1/*okay*/ : 0/*failure*/;
    } else {
        return s_OctalToNum  (val, ptr, len) ?  1/*okay*/ : 0/*failure*/;
    }
}


static mode_t s_TarToMode(TTarMode value)
{
    mode_t mode = (
#ifdef S_ISUID
                   (value & fTarSetUID   ? S_ISUID  : 0) |
#endif
#ifdef S_ISGID
                   (value & fTarSetGID   ? S_ISGID  : 0) |
#endif
#ifdef S_ISVTX
                   (value & fTarSticky   ? S_ISVTX  : 0) |
#endif
#if   defined(S_IRUSR)
                   (value & fTarURead    ? S_IRUSR  : 0) |
#elif defined(S_IREAD)
                   (value & fTarURead    ? S_IREAD  : 0) |
#endif
#if   defined(S_IWUSR)
                   (value & fTarUWrite   ? S_IWUSR  : 0) |
#elif defined(S_IWRITE)
                   (value & fTarUWrite   ? S_IWRITE : 0) |
#endif
#if   defined(S_IWUSR)
                   (value & fTarUExecute ? S_IXUSR  : 0) |
#elif defined(S_IEXEC)
                   (value & fTarUExecute ? S_IEXEC  : 0) |
#endif
#ifdef S_IRGRP
                   (value & fTarGRead    ? S_IRGRP  : 0) |
#endif
#ifdef S_IWGRP
                   (value & fTarGWrite   ? S_IWGRP  : 0) |
#endif
#ifdef S_IXGRP
                   (value & fTarGExecute ? S_IXGRP  : 0) |
#endif
#ifdef S_IROTH
                   (value & fTarORead    ? S_IROTH  : 0) |
#endif
#ifdef S_IWOTH
                   (value & fTarOWrite   ? S_IWOTH  : 0) |
#endif
#ifdef S_IXOTH
                   (value & fTarOExecute ? S_IXOTH  : 0) |
#endif
                   0);
    return mode;
}


static TTarMode s_ModeToTar(mode_t mode)
{
    // Keep in mind that the mode may be extracted on a different platform
    TTarMode value = (
#ifdef S_ISUID
                      (mode & S_ISUID  ? fTarSetUID   : 0) |
#endif
#ifdef S_ISGID
                      (mode & S_ISGID  ? fTarSetGID   : 0) |
#endif
#ifdef S_ISVTX
                      (mode & S_ISVTX  ? fTarSticky   : 0) |
#endif
#if   defined(S_IRUSR)
                      (mode & S_IRUSR  ? fTarURead    : 0) |
#elif defined(S_IREAD)
                      (mode & S_IREAD  ? fTarURead    : 0) |
#endif
#if   defined(S_IWUSR)
                      (mode & S_IWUSR  ? fTarUWrite   : 0) |
#elif defined(S_IWRITE)
                      (mode & S_IWRITE ? fTarUWrite   : 0) |
#endif
#if   defined(S_IXUSR)
                      (mode & S_IXUSR  ? fTarUExecute : 0) |
#elif defined(S_IEXEC)
                      (mode & S_IEXEC  ? fTarUExecute : 0) |
#endif
#if   defined(S_IRGRP)
                      (mode & S_IRGRP  ? fTarGRead    : 0) |
#elif defined(S_IREAD)
                      // emulate read permission when file is readable
                      (mode & S_IREAD  ? fTarGRead    : 0) |
#endif
#ifdef S_IWGRP
                      (mode & S_IWGRP  ? fTarGWrite   : 0) |
#endif
#ifdef S_IXGRP
                      (mode & S_IXGRP  ? fTarGExecute : 0) |
#endif
#if   defined(S_IROTH)
                      (mode & S_IROTH  ? fTarORead    : 0) |
#elif defined(S_IREAD)
                      // emulate read permission when file is readable
                      (mode & S_IREAD  ? fTarORead    : 0) |
#endif
#ifdef S_IWOTH
                      (mode & S_IWOTH  ? fTarOWrite   : 0) |
#endif
#ifdef S_IXOTH
                      (mode & S_IXOTH  ? fTarOExecute : 0) |
#endif
                      0);
#if defined(S_IFMT)  ||  defined(_S_IFMT)
    TTarMode mask = (TTarMode) mode;
#  ifdef S_IFMT
    mask &=  S_IFMT;
#  else
    mask &= _S_IFMT;
#  endif
    if (!(mask & 07777)) {
        value |= mask;
    }
#endif
    return value;
}


static size_t s_Length(const char* ptr, size_t maxsize)
{
    const char* pos = (const char*) memchr(ptr, '\0', maxsize);
    return pos ? (size_t)(pos - ptr) : maxsize;
}


//////////////////////////////////////////////////////////////////////////////
//
// Constants / macros / typedefs
//

/// Round up to the nearest multiple of kBlockSize:
//#define ALIGN_SIZE(size)   SIZE_OF(BLOCK_OF(size + (BLOCK_SIZE-1)))
#define ALIGN_SIZE(size) (((size) + (BLOCK_SIZE-1)) & ~(BLOCK_SIZE-1))
#define OFFSET_OF(size)  ( (size)                   &  (BLOCK_SIZE-1))
#define BLOCK_OF(pos)    ((pos) >> 9)
#define SIZE_OF(blk)     ((blk) << 9)

/// Tar block size (512 bytes)
/// NB: Using Uint8 here is bad because we want to have "regular sizes" as
/// size_t.  Using size_t here is bad because (see ALIGN_SIZE above) it
/// can cause 32-bit truncation of the result (as in ~(BLOCK_SIZE-1) being
/// zero-filled in the higher-order bits if sizeof(size_t)==4) in expressions
/// where mixed with wider Uint8's.  So the macro is the best bet!
#define BLOCK_SIZE       SIZE_OF(1)


/// Recognized TAR formats
enum ETar_Format {
    eTar_Unknown = 0,
    eTar_Legacy  = 1,
    eTar_OldGNU  = 2,
    eTar_Ustar   = 4,
    eTar_Posix   = 5
};


/// POSIX "ustar" tar archive member header
typedef struct SHeader {     // byte offset
    char name[100];          //   0
    char mode[8];            // 100
    char uid[8];             // 108
    char gid[8];             // 116
    char size[12];           // 124
    char mtime[12];          // 136
    char checksum[8];        // 148
    char typeflag[1];        // 156
    char linkname[100];      // 157
    char magic[6];           // 257
    char version[2];         // 263
    char uname[32];          // 265
    char gname[32];          // 297
    char devmajor[8];        // 329
    char devminor[8];        // 337
    union {                  // 345
        char prefix[155];    // NB: not valid with old GNU format (no need)
        struct {             // NB:                old GNU format only
            char atime[12];
            char ctime[12];  // 357
        } gt;
    };                       // 500
    // NCBI in last 4 bytes  // 508
} SHeader;


/// Block as a header.
union TBlock {
    char    buffer[BLOCK_SIZE];
    SHeader header;
};


static bool s_TarChecksum(TBlock* block, bool isgnu)
{
    SHeader* h = &block->header;
    size_t len = sizeof(h->checksum) - (isgnu ? 2 : 1);

    // Compute the checksum
    memset(h->checksum, ' ', sizeof(h->checksum));
    unsigned long checksum = 0;
    const unsigned char* p = (const unsigned char*) block->buffer;
    for (size_t i = 0; i < sizeof(block->buffer); i++) {
        checksum += *p++;
    }
    // ustar:       '\0'-terminated checksum
    // GNU special: 6 digits, then '\0', then a space [already in place]
    if (!s_NumToOctal(checksum, h->checksum, len)) {
        return false;
    }
    h->checksum[len] = '\0';
    return true;
}



//////////////////////////////////////////////////////////////////////////////
//
// CTarEntryInfo
//

TTarMode CTarEntryInfo::GetMode(void) const
{
    // Raw tar mode gets returned here (as kept in the info)
    return (TTarMode)(m_Stat.st_mode & 07777);
}


void CTarEntryInfo::GetMode(CDirEntry::TMode*            usr_mode,
                            CDirEntry::TMode*            grp_mode,
                            CDirEntry::TMode*            oth_mode,
                            CDirEntry::TSpecialModeBits* special_bits) const
{
    TTarMode mode = GetMode();

    // User
    if (usr_mode) {
        *usr_mode = ((mode & fTarURead    ? CDirEntry::fRead    : 0) |
                     (mode & fTarUWrite   ? CDirEntry::fWrite   : 0) |
                     (mode & fTarUExecute ? CDirEntry::fExecute : 0));
    }

    // Group
    if (grp_mode) {
        *grp_mode = ((mode & fTarGRead    ? CDirEntry::fRead    : 0) |
                     (mode & fTarGWrite   ? CDirEntry::fWrite   : 0) |
                     (mode & fTarGExecute ? CDirEntry::fExecute : 0));
    }

    // Others
    if (oth_mode) {
        *oth_mode = ((mode & fTarORead    ? CDirEntry::fRead    : 0) |
                     (mode & fTarOWrite   ? CDirEntry::fWrite   : 0) |
                     (mode & fTarOExecute ? CDirEntry::fExecute : 0));
    }

    // Special bits
    if (special_bits) {
        *special_bits = ((mode & fTarSetUID ? CDirEntry::fSetUID : 0) |
                         (mode & fTarSetGID ? CDirEntry::fSetGID : 0) |
                         (mode & fTarSticky ? CDirEntry::fSticky : 0));
    }

    return;
}


unsigned int CTarEntryInfo::GetMajor(void) const
{
#ifdef major
    if (m_Type == eCharDev  ||  m_Type == eBlockDev) {
        return major(m_Stat.st_rdev);
    }
#else
    if (sizeof(int) >= 4  &&  sizeof(m_Stat.st_rdev) >= 4) {
        return (*((unsigned int*) &m_Stat.st_rdev) >> 16) & 0xFFFF;
    }
#endif // major
    return (unsigned int)(-1);
}


unsigned int CTarEntryInfo::GetMinor(void) const
{
#ifdef minor
    if (m_Type == eCharDev  ||  m_Type == eBlockDev) {
        return minor(m_Stat.st_rdev);
    }
#else
    if (sizeof(int) >= 4  &&  sizeof(m_Stat.st_rdev) >= 4) {
        return *((unsigned int*) &m_Stat.st_rdev) & 0xFFFF;
    }
#endif // minor
    return (unsigned int)(-1);
}


bool CTarEntryInfo::operator==(const CTarEntryInfo& info) const
{
    return (m_Type       == info.m_Type                        &&
            m_Name       == info.m_Name                        &&
            m_LinkName   == info.m_LinkName                    &&
            m_UserName   == info.m_UserName                    &&
            m_GroupName  == info.m_GroupName                   &&
            m_HeaderSize == info.m_HeaderSize                  &&
            memcmp(&m_Stat,&info.m_Stat, sizeof(m_Stat)) == 0  &&
            m_Pos        == info.m_Pos ? true : false);
}


static string s_ModeAsString(TTarMode mode)
{
    char  buf[9];
    memset(buf, '-', sizeof(buf));

    char* usr = buf;
    char* grp = usr + 3;
    char* oth = grp + 3;

    if (mode & fTarURead) {
        usr[0] = 'r';
    }
    if (mode & fTarUWrite) {
        usr[1] = 'w';
    }
    if (mode & fTarUExecute) {
        usr[2] = mode & fTarSetUID ? 's' : 'x';
    } else if (mode & fTarSetUID) {
        usr[2] = 'S';
    }
    if (mode & fTarGRead) {
        grp[0] = 'r';
    }
    if (mode & fTarGWrite) {
        grp[1] = 'w';
    }
    if (mode & fTarGExecute) {
        grp[2] = mode & fTarSetGID ? 's' : 'x';
    } else if (mode & fTarSetGID) {
        grp[2] = 'S';
    }
    if (mode & fTarORead) {
        oth[0] = 'r';
    }
    if (mode & fTarOWrite) {
        oth[1] = 'w';
    }
    if (mode & fTarOExecute) {
        oth[2] = mode & fTarSticky ? 't' : 'x';
    } else if (mode & fTarSticky) {
        oth[2] = 'T';
    }

    return string(buf, sizeof(buf));
}


static char s_TypeAsChar(CTarEntryInfo::EType type)
{
    switch (type) {
    case CTarEntryInfo::eFile:
    case CTarEntryInfo::eHardLink:
        return '-';
    case CTarEntryInfo::eSymLink:
        return 'l';
    case CTarEntryInfo::eDir:
        return 'd';
    case CTarEntryInfo::ePipe:
        return 'p';
    case CTarEntryInfo::eCharDev:
        return 'c';
    case CTarEntryInfo::eBlockDev:
        return 'b';
    default:
        break;
    }
    return '?';
}


static string s_UserGroupAsString(const CTarEntryInfo& info)
{
    string user(info.GetUserName());
    if (user.empty()) {
        NStr::UIntToString(user, info.GetUserId());
    }
    string group(info.GetGroupName());
    if (group.empty()) {
        NStr::UIntToString(group, info.GetGroupId());
    }
    return user + "/" + group;
}


static string s_MajorMinor(unsigned int n)
{
    return n == (unsigned int)(-1) ? string("?") : NStr::UIntToString(n);
}


static string s_SizeOrMajorMinor(const CTarEntryInfo& info)
{
    if (info.GetType() == CTarEntryInfo::eCharDev  ||
        info.GetType() == CTarEntryInfo::eBlockDev) {
        unsigned int major = info.GetMajor();
        unsigned int minor = info.GetMinor();
        return s_MajorMinor(major) + ", " + s_MajorMinor(minor);
    }
    return NStr::UInt8ToString(info.GetSize());
}


ostream& operator << (ostream& os, const CTarEntryInfo& info)
{
    CTime mtime(info.GetModificationTime());
    os << s_TypeAsChar(info.GetType())
       << s_ModeAsString(info.GetMode())        << ' '
       << setw(17) << s_UserGroupAsString(info) << ' '
       << setw(10) << s_SizeOrMajorMinor(info)  << ' '
       << mtime.ToLocalTime().AsString(" Y-M-D h:m:s ")
       << info.GetName();
    if (info.GetType() == CTarEntryInfo::eSymLink  ||
        info.GetType() == CTarEntryInfo::eHardLink) {
        os << " -> " << info.GetLinkName();
    }
    return os;
}



//////////////////////////////////////////////////////////////////////////////
//
// Debugging utilities
//

static string s_OSReason(int x_errno)
{
    const char* strerr = x_errno ? strerror(x_errno) : 0;
    return strerr  &&  *strerr ? string(": ") + strerr : kEmptyStr;
}


static string s_PositionAsString(const string& file, Uint8 pos, size_t recsize,
                                 const string& entryname)
{
    _ASSERT(!OFFSET_OF(pos));
    _ASSERT(!OFFSET_OF(recsize));
    _ASSERT(recsize >= BLOCK_SIZE);
    string result;
    if (!file.empty()) {
        CDirEntry temp(file);
        result = temp.GetName() + ": ";
    }
    result += "At record " + NStr::UInt8ToString(pos / recsize);
    if (recsize != BLOCK_SIZE) {
        result +=
            ", block " + NStr::UInt8ToString(BLOCK_OF(pos % recsize)) +
            " [thru #" + NStr::UInt8ToString(BLOCK_OF(pos),
                                             NStr::fWithCommas) + ']';
    }
    if (!entryname.empty()) {
        result += ", while in '" + entryname + '\'';
    }
    return result + ":\n";
}


static string s_AddressAsString(size_t offset)
{
    // All addresses but the base are greater than 100:
    // do the quick and dirty formatting here
    return offset ? NStr::UIntToString((unsigned long) offset) : "000";
}


static bool memcchr(const char* s, char c, size_t len)
{
    for (size_t i = 0;  i < len;  i++) {
        if (s[i] != c)
            return true;
    }
    return false;
}


static string s_Printable(const char* field, size_t maxsize, bool text)
{
    if (!text  &&  maxsize > 1  &&  !*field) {
        field++, maxsize--;
    }
    size_t len = s_Length(field, maxsize);
    return NStr::PrintableString(string(field,
                                        memcchr(field + len, '\0',
                                                maxsize - len)
                                        ? maxsize
                                        : len));
}


#if !defined(__GNUC__)  &&  !defined(offsetof)
#  define offsetof(T, F) ((char*) &(((T*) 0)->F) - (char*) 0)
#endif


#define _STR(s) #s

#define TAR_PRINTABLE(field, text)                                      \
    "@" + s_AddressAsString((size_t) offsetof(SHeader, field)) +        \
    '[' + _STR(field) + "]:" + string(10 - sizeof(_STR(field)), ' ') +  \
    '"' + s_Printable(h->field, sizeof(h->field), text  ||  ex) + '"'

static string s_DumpHeader(const SHeader* h, ETar_Format fmt, bool ex = false)
{
    string dump;
    Uint8 val;
    int ok;

    dump += TAR_PRINTABLE(name, true) + '\n';

    ok = s_OctalToNum(val, h->mode, sizeof(h->mode));
    dump += TAR_PRINTABLE(mode, !ok);
    if (ok  &&  val) {
        dump += " [" + s_ModeAsString((TTarMode) val) + ']';
    }
    dump += '\n';

    ok = s_DecodeUint8(val, h->uid, sizeof(h->uid));
    dump += TAR_PRINTABLE(uid, ok <= 0);
    if (ok  &&  (ok < 0  || val > 7)) {
        dump += " [" + NStr::UInt8ToString(val) + ']';
        if (ok < 0) {
            dump += " (base-256)";
        }
    }
    dump += '\n';
    
    ok = s_DecodeUint8(val, h->gid, sizeof(h->gid));
    dump += TAR_PRINTABLE(gid, ok <= 0);
    if (ok  &&  (ok < 0  ||  val > 7)) {
        dump += " [" + NStr::UInt8ToString(val) + ']';
        if (ok < 0) {
            dump += " (base-256)";
        }
    }
    dump += '\n';

    ok = s_DecodeUint8(val, h->size, sizeof(h->size));
    dump += TAR_PRINTABLE(size, ok <= 0);
    if (ok  &&  (ok < 0  ||  val > 7)) {
        dump += " [" + NStr::UInt8ToString(val) + ']';
        if (ok < 0) {
            dump += " (base-256)";
        }
    }
    dump += '\n';

    ok = s_OctalToNum(val, h->mtime, sizeof(h->mtime));
    dump += TAR_PRINTABLE(mtime, !ok);
    if (ok  &&  val) {
        CTime mtime((time_t) val);
        ok = (Uint8) mtime.GetTimeT() == val ? true : false;
        if (ok  ||  val > 7) {
            dump += (" ["
                     + (val > 7 ? NStr::UInt8ToString(val) + ", "        : "")
                     + (ok ? mtime.ToLocalTime().AsString("Y-M-D h:m:s") : "")
                     + ']');
        }
    }
    dump += '\n';

    ok = s_OctalToNum(val, h->checksum, sizeof(h->checksum));
    dump += TAR_PRINTABLE(checksum, !ok) + '\n';

    // Classify to the extent possible to help debug the problem (if any)
    dump += TAR_PRINTABLE(typeflag, true);
    ok = false;
    const char* tname = 0;
    switch (h->typeflag[0]) {
    case '\0':
    case '0':
        ok = true;
        if (!(fmt & eTar_Ustar)  &&  fmt != eTar_OldGNU) {
            size_t namelen = s_Length(h->name, sizeof(h->name));
            if (namelen  &&  h->name[namelen - 1] == '/') {
                tname = "legacy regular entry (dir)" + (h->typeflag[0]? 7 : 0);
                break;
            }
        }
        tname = "legacy regular entry (file)" + (h->typeflag[0] ? 7 : 0);
        break;
    case '1':
        ok = true;
#ifdef NCBI_OS_UNIX
        tname = "hard link";
#else
        tname = "hard link - not FULLY supported";
#endif // NCBI_OS_UNIX
        break;
    case '2':
        ok = true;
#ifdef NCBI_OS_UNIX
        tname = "symbolic link";
#else
        tname = "symbolic link - not FULLY supported";
#endif // NCBI_OS_UNIX
        break;
    case '3':
#ifdef NCBI_OS_UNIX
        ok = true;
#endif // NCBI_OS_UNIX
        tname = "character device";
        break;
    case '4':
#ifdef NCBI_OS_UNIX
        ok = true;
#endif // NCBI_OS_UNIX
        tname = "block device";
        break;
    case '5':
        ok = true;
        tname = "directory";
        break;
    case '6':
#ifdef NCBI_OS_UNIX
        ok = true;
#endif // NCBI_OS_UNIX
        tname = "FIFO";
        break;
    case '7':
        tname = "contiguous";
        break;
    case 'g':
        tname = "global extended header";
        break;
    case 'x':
    case 'X':
        if (fmt & eTar_Ustar) {
            ok = true;
            if (h->typeflag[0] == 'x') {
                tname = "extended (POSIX 1003.1-2001 [PAX]) header"
                    " - not FULLY supported";
            } else {
                tname = "extended (POSIX 1003.1-2001 [PAX] by Sun) header"
                    " - not FULLY supported";
            }
        } else {
            tname = "extended header";
        }
        break;
    case 'A':
        tname = "Solaris ACL";
        break;
    case 'D':
        if (fmt == eTar_OldGNU) {
            tname = "GNU extension: directory dump";
        }
        break;
    case 'E':
        tname = "Solaris extended attribute file";
        break;
    case 'I':
        tname = "inode only";
        break;
    case 'K':
        if (fmt == eTar_OldGNU) {
            ok = true;
            tname = "GNU extension: long link";
        }
        break;
    case 'L':
        if (fmt == eTar_OldGNU) {
            ok = true;
            tname = "GNU extension: long name";            
        }
        break;
    case 'M':
        if (fmt == eTar_OldGNU) {
            tname = "GNU extension: multivolume archive";
        }
        break;
    case 'N':
        if (fmt == eTar_OldGNU) {
            tname = "GNU extension: long filename";
        }
        break;
    case 'S':
        if (fmt == eTar_OldGNU) {
            tname = "GNU extension: sparse file";
        }
        break;
    case 'V':
        if (fmt == eTar_OldGNU) {
            tname = "GNU extension: volume header";
        }
        break;
    default:
        break;
    }
    if (!tname  &&  h->typeflag[0] >= 'A'  &&  h->typeflag[0] <= 'Z') {
        tname = "user-defined extension";
    }
    dump += (" [" + string(tname ? tname : "reserved") +
             (ok
              ? "]\n"
              : " -- NOT SUPPORTED]\n"));

    dump += TAR_PRINTABLE(linkname, true) + '\n';

    switch (fmt) {
    case eTar_Legacy:  // NCBI does not ever write this header
        tname = "legacy (V7)";
        break;
    case eTar_OldGNU:
        if (!NStr::strncasecmp((const char*) h + BLOCK_SIZE - 4, "NCBI", 4)) {
            tname = "old GNU (NCBI)";
        } else {
            tname = "old GNU";
        }
        break;
    case eTar_Ustar:
        if (!NStr::strncasecmp((const char*) h + BLOCK_SIZE - 4, "NCBI", 4)) {
            tname = "ustar (NCBI)";
        } else {
            tname = "ustar";
        }
        break;
    case eTar_Posix:  // aka "pax"
        if (!NStr::strncasecmp((const char*) h + BLOCK_SIZE - 4, "NCBI", 4)) {
            tname = "posix (NCBI)";
        } else {
            tname = "posix";
        }
        break;
    default:
        tname = 0;
        break;
    }
    dump += TAR_PRINTABLE(magic, true);
    if (tname) {
        dump += " [" + string(tname) + ']';
    }
    dump += '\n';

    dump += TAR_PRINTABLE(version, true);

    if (fmt != eTar_Legacy) {
        dump += '\n';

        dump += TAR_PRINTABLE(uname, true) + '\n';

        dump += TAR_PRINTABLE(gname, true) + '\n';

        ok = s_OctalToNum(val, h->devmajor, sizeof(h->devmajor));
        dump += TAR_PRINTABLE(devmajor, !ok);
        if (ok  &&  val > 7) {
            dump += " [" + NStr::UInt8ToString(val) + ']';
        }
        dump += '\n';

        ok = s_OctalToNum(val, h->devminor, sizeof(h->devminor));
        dump += TAR_PRINTABLE(devminor, !ok);
        if (ok  &&  val > 7) {
            dump += " [" + NStr::UInt8ToString(val) + ']';
        }
        dump += '\n';

        if (fmt == eTar_OldGNU) {
            ok = s_OctalToNum(val, h->gt.atime, sizeof(h->gt.atime));
            dump += TAR_PRINTABLE(gt.atime, !ok);
            if (ok  &&  val) {
                CTime atime((time_t) val);
                ok = (Uint8) atime.GetTimeT() == val ? true : false;
                if (ok  ||  val > 7) {
                    dump += (" ["
                             + (val > 7 ? NStr::UInt8ToString(val) + ", " : "")
                             + (ok
                                ? atime.ToLocalTime().AsString("Y-M-D h:m:s")
                                : "")
                             + ']');
                }
            }
            dump += '\n';

            ok = s_OctalToNum(val, h->gt.ctime, sizeof(h->gt.ctime));
            dump += TAR_PRINTABLE(gt.ctime, !ok);
            if (ok  &&  val) {
                CTime ctime((time_t) val);
                ok = (Uint8) ctime.GetTimeT() == val ? true : false;
                if (ok  ||  val > 7) {
                    dump += (" ["
                             + (val > 7 ? NStr::UInt8ToString(val) + ", " : "")
                             + (ok
                                ? ctime.ToLocalTime().AsString("Y-M-D h:m:s")
                                : "")
                             + ']');
                }
            }
            tname = h->gt.ctime + sizeof(h->gt.ctime);
        } else {
            dump += TAR_PRINTABLE(prefix, true);
            tname = h->prefix + sizeof(h->prefix);
        }
    } else {
        tname = h->version + sizeof(h->version);
    }

    size_t n = 0;
    while (&tname[n] < (const char*) h + BLOCK_SIZE) {
        if (tname[n]) {
            size_t offset = (size_t)(&tname[n] - (const char*) h);
            size_t len = BLOCK_SIZE - offset;
            if (len & ~0xF) {  // len > 16
                len = 0x10;    // len = 16
            }
            const char* e = (const char*) memchr(&tname[n], '\0', len);
            if (e) {
                len = (size_t)(e - &tname[n]);
                ok = s_DecodeUint8(val, &tname[n], len);
            } else {
                if (len  > (offset & 0xF)) {
                    len -= (offset & 0xF);
                }
                ok = false;
            }
            _ASSERT(len);
            dump += "\n@" + s_AddressAsString(offset) + ':' + string(11, ' ')
                + '"' + NStr::PrintableString(string(&tname[n], len)) + '"';
            if (ok) {
                CTime time((time_t) val);
                bool okaytime = (Uint8) time.GetTimeT() == val;
                if (ok < 0  ||  val > 7  ||  okaytime) {
                    dump += " [";
                    if (ok < 0  ||  val > 7) {
                        dump += NStr::UInt8ToString(val);
                    }
                    if (ok < 0) {
                        dump += "] (base-256)";
                    } else if (okaytime) {
                        if (val > 7) {
                            dump += ", ";
                        }
                        dump += time.ToLocalTime().AsString("Y-M-D h:m:s]");
                    } else {
                        dump += ']';
                    }
                }
            }
            n += len;
        } else {
            n++;
        }
    }

    return dump;
}

#undef TAR_PRINTABLE

#undef _STR


//////////////////////////////////////////////////////////////////////////////
//
// CTar
//

CTar::CTar(const string& filename, size_t blocking_factor)
    : m_FileName(filename),
      m_FileStream(new CNcbiFstream),
      m_OpenMode(eNone),
      m_Stream(0),
      m_BufferSize(SIZE_OF(blocking_factor)),
      m_BufferPos(0),
      m_StreamPos(0),
      m_BufPtr(0),
      m_Buffer(0),
      m_Mask(0),
      m_MaskOwned(eNoOwnership),
      m_Modified(false),
      m_Bad(false),
      m_Flags(fDefault)
{
    x_Init();
}


CTar::CTar(CNcbiIos& stream, size_t blocking_factor)
    : m_FileName(kEmptyStr),
      m_FileStream(0),
      m_OpenMode(eNone),
      m_Stream(&stream),
      m_BufferSize(SIZE_OF(blocking_factor)),
      m_BufferPos(0),
      m_StreamPos(0),
      m_BufPtr(0),
      m_Buffer(0),
      m_Mask(0),
      m_MaskOwned(eNoOwnership),
      m_Modified(false),
      m_Bad(false),
      m_Flags(fDefault)
{
    x_Init();
}


CTar::~CTar()
{
    // Close stream(s)
    x_Flush(true/*no_throw*/);
    x_Close();
    delete m_FileStream;

    // Delete owned file name masks
    UnsetMask();

    // Delete buffer
    delete[] m_BufPtr;
}


#define TAR_THROW(who, errcode, message)                                \
    NCBI_THROW(CTarException, errcode,                                  \
               s_PositionAsString(who->m_FileName, who->m_StreamPos,    \
                                  who->m_BufferSize,                    \
                                  who->m_Current.GetName()) + (message))

#define TAR_THROW_EX(who, errcode, message, h, fmt)                     \
    TAR_THROW(who, errcode,                                             \
              who->m_Flags & fDumpBlockHeaders                          \
              ? string(message) + ":\n" + s_DumpHeader(h, fmt, true)    \
              : string(message))

#define TAR_POST(subcode, severity, message)                            \
    ERR_POST_X(subcode, severity <<                                     \
               s_PositionAsString(m_FileName, m_StreamPos, m_BufferSize,\
                                  m_Current.GetName()) + (message))


void CTar::x_Init(void)
{
    size_t pagesize = (size_t) GetVirtualMemoryPageSize();
    if (!pagesize) {
        pagesize = 4096; // reasonable default
    }
    size_t pagemask = pagesize - 1;
    // Assume that the page size is a power of 2
    _ASSERT((pagesize & pagemask) == 0);
    m_BufPtr = new char[m_BufferSize + pagemask];
    // Make m_Buffer page-aligned
    m_Buffer = m_BufPtr +
        ((((size_t) m_BufPtr + pagemask) & ~pagemask) - (size_t) m_BufPtr);
}


void CTar::x_Flush(bool no_throw)
{
    m_Current.m_Name.erase();
    if (!m_Stream  ||  !m_OpenMode  ||  !m_Modified) {
        return;
    }
    _ASSERT(m_OpenMode == eRW);
    _ASSERT(m_BufferSize >= m_BufferPos);
    size_t pad = m_BufferSize - m_BufferPos;
    // Assure proper blocking factor and pad the archive as necessary
    memset(m_Buffer + m_BufferPos, 0, pad);
    x_WriteArchive(pad, no_throw ? (const char*)(-1L) : 0);
    _ASSERT(!(m_BufferPos % m_BufferSize)); // m_BufferSize if write error
    _ASSERT(!m_Bad == !m_BufferPos);
    if (!m_Bad
        &&  (pad < BLOCK_SIZE || pad - OFFSET_OF(pad) < (BLOCK_SIZE << 1))) {
        // Write EOT (two zero blocks), if have not already done so by padding
        memset(m_Buffer, 0, m_BufferSize - pad);
        x_WriteArchive(m_BufferSize, no_throw ? (const char*)(-1L) : 0);
        _ASSERT(!(m_BufferPos % m_BufferSize)  &&  !m_Bad == !m_BufferPos);
        if (!m_Bad  &&  m_BufferSize == BLOCK_SIZE) {
            x_WriteArchive(BLOCK_SIZE, no_throw ? (const char*)(-1L) : 0);
            _ASSERT(!(m_BufferPos % m_BufferSize)  &&  !m_Bad == !m_BufferPos);
        }
    }
    if (!m_Bad  &&  m_Stream->rdbuf()->PUBSYNC() != 0) {
        int x_errno = errno;
        if (!no_throw) {
            TAR_THROW(this, eWrite,
                      "Archive flush failed" + s_OSReason(x_errno));
        }
        TAR_POST(83, Error,
                 "Archive flush failed" + s_OSReason(x_errno));
    }
    m_Modified = false;
}


void CTar::x_Close(void)
{
    if (!m_OpenMode) {
        return;
    }
    if (m_FileStream  &&  m_FileStream->is_open()) {
        m_FileStream->close();
        m_Stream = 0;
        m_Bad = false;
    }
    m_OpenMode  = eNone;
    m_BufferPos = 0;
    m_StreamPos = 0;
}


auto_ptr<CTar::TEntries> CTar::x_Open(EAction action)
{
    _ASSERT(action);
    // We can only open a named file here, and if an external stream
    // is being used as an archive, it must be explicitly repositioned by
    // user's code (outside of this class) before each archive operation.
    if (!m_FileStream) {
        if (m_Modified  &&  action != eAppend) {
            TAR_POST(1, Warning,
                     "Pending changes may be discarded"
                     " upon reopen of in-stream archive");
            m_Modified = false;
        }
        if (action != eInternal) {
            m_BufferPos = 0;
            m_StreamPos = 0;
        }
        m_Current.m_Name.erase();
        if (m_Bad || !m_Stream || !m_Stream->good() || !m_Stream->rdbuf()) {
            m_OpenMode = eNone;
            TAR_THROW(this, eOpen,
                      "Bad IO stream provided for archive");
        } else {
            m_OpenMode = EOpenMode(int(action) & eRW);
        }
#ifdef NCBI_OS_MSWIN
        if (m_Stream == &cin) {
            HANDLE handle = (HANDLE) _get_osfhandle(_fileno(stdin));
            if (GetFileType(handle) != FILE_TYPE_DISK) {
                m_Flags |= fSlowSkipWithRead;
            }
        }
#endif // NCBI_OS_MSWIN
    } else {
        EOpenMode mode = EOpenMode(int(action) & eRW);
        _ASSERT(mode != eNone);
        if (mode != eWO  &&  action != eAppend) {
            x_Flush();
        } else {
            m_Current.m_Name.erase();
        }
        if (mode == eWO  ||  m_OpenMode < mode) {
            x_Close();
            switch (mode) {
            case eWO:
                // WO access
                m_FileStream->open(m_FileName.c_str(),
                                   IOS_BASE::out    |
                                   IOS_BASE::binary | IOS_BASE::trunc);
                break;
            case eRO:
                // RO access
                m_FileStream->open(m_FileName.c_str(),
                                   IOS_BASE::in     |
                                   IOS_BASE::binary);
                break;
            case eRW:
                // RW access
                m_FileStream->open(m_FileName.c_str(),
                                   IOS_BASE::in     | IOS_BASE::out |
                                   IOS_BASE::binary);
                break;
            default:
                _TROUBLE;
                break;
            }
            if (!m_FileStream->good()) {
                int x_errno = errno;
                TAR_THROW(this, eOpen,
                          "Cannot open archive '" + m_FileName + '\''
                          + s_OSReason(x_errno));
            }
            m_Stream = m_FileStream;
            m_Bad = false;
        }

        if (m_OpenMode) {
            _ASSERT(action != eCreate);
            if (action != eAppend  &&  action != eInternal) {
                m_BufferPos = 0;
                m_StreamPos = 0;
                m_FileStream->seekg(0, IOS_BASE::beg);
            }
        } else {
            if (action == eAppend  &&  !m_Modified) {
                // There may be an extra and unnecessary archive scanning
                // if Append() follows Update() that caused no modifications;
                // but there is no way to distinguish this, currently :-/
                // Also, this sequence should be a real rarity in practice.
                x_ReadAndProcess(eAppend);  // positions at logical EOF
            }
            m_OpenMode = mode;
        }
    }
    _ASSERT(m_Stream);
    _ASSERT(m_Stream->rdbuf());

    if (int(action) & ((eList | eExtract | eInternal) & ~eRW)) {
        return x_ReadAndProcess(action);
    } else {
        return auto_ptr<TEntries>(0);
    }
}


auto_ptr<CTar::TEntries> CTar::Extract(void)
{
    auto_ptr<TEntries> done = x_Open(eExtract);

    // Restore attributes of "postponed" directory entries
    if (m_Flags & fPreserveAll) {
        ITERATE(TEntries, i, *done.get()) {
            if (i->GetType() == CTarEntryInfo::eDir) {
                x_RestoreAttrs(*i);
            }
        }
    }

    return done;
}


const CTarEntryInfo* CTar::GetNextEntryInfo(void)
{
    if (m_OpenMode) {
        x_Skip(BLOCK_OF(m_Current.GetPosition(CTarEntryInfo::ePos_Data)
                        + ALIGN_SIZE(m_Current.GetSize()) - m_StreamPos));
    }

    auto_ptr<TEntries> temp = x_Open(eInternal);
    _ASSERT(temp.get()  &&  temp->size() < 2);
    if (temp->size() < 1) {
        return 0;
    }

    _ASSERT(m_Current == *temp->begin());
    return &m_Current;
}


Uint8 CTar::EstimateArchiveSize(const TFiles& files) const
{
    Uint8 result = 0;

    ITERATE(TFiles, it, files) {
        // Count in the file size
        result += BLOCK_SIZE/*header*/ + ALIGN_SIZE(it->second);

        // Count in the long name (if any)
        string path = x_ToFilesystemPath(it->first);
        string name = x_ToArchiveName(path);
        size_t namelen = name.length() + 1;
        if (namelen > sizeof(((SHeader*) 0)->name)) {
            result += BLOCK_SIZE/*long name header*/ + ALIGN_SIZE(namelen);
        }
    }
    if (result) {
        result += BLOCK_SIZE << 1; // EOT
        Uint8 incomplete = result % m_BufferSize;
        if (incomplete) {
            result += m_BufferSize - incomplete;
        }
    }

    return result;
}


void CTar::SetBaseDir(const string& dirname)
{
    m_BaseDir = CDirEntry::AddTrailingPathSeparator(dirname);
#ifdef NCBI_OS_MSWIN
    // Replace backslashes with forward slashes
    NStr::ReplaceInPlace(m_BaseDir, "\\", "/");
#endif // NCBI_OS_MSWIN
}


// Return a pointer to buffer, which is always block-aligned,
// and reflect the number of bytes available via the parameter.
const char* CTar::x_ReadArchive(size_t& n)
{
    _ASSERT(!OFFSET_OF(m_BufferPos));
    _ASSERT(n != 0);
    size_t nread;
    if (!m_BufferPos) {
        nread = 0;
        do {
#ifdef NCBI_COMPILER_MIPSPRO
            // Work around a bug in MIPSPro 7.3's streambuf::xsgetn()
            istream* is = dynamic_cast<istream*>(m_Stream);
            _ASSERT(is);
            is->read(             m_Buffer     + nread,
                     (streamsize)(m_BufferSize - nread));
            streamsize xread = is->gcount();
            is->clear();
#else
            streamsize xread = m_Stream->rdbuf()
                ->sgetn(             m_Buffer     + nread,
                        (streamsize)(m_BufferSize - nread));
#endif // NCBI_COMPILER_MIPSPRO
            if (xread <= 0) {
                break;
            }
            nread += xread;
        } while (nread < m_BufferSize);
        size_t gap = m_BufferSize - nread;
        if (gap < m_BufferSize) {
            memset(m_Buffer + nread, 0, gap);
        } else {
            return 0/*EOF*/;
        }
    } else {
        nread = m_BufferSize - m_BufferPos;
    }
    if (n > nread) {
        n = nread;
    }
    size_t xpos = m_BufferPos;
    m_BufferPos += ALIGN_SIZE(n);
    m_BufferPos %= m_BufferSize;
    return m_Buffer + xpos;
}


// All partial internal (i.e. in-buffer) block writes are _not_ block-aligned.
void CTar::x_WriteArchive(size_t nwrite, const char* src)
{
    if (!nwrite  ||  m_Bad) {
        return;
    }
    m_Modified = true;
    do {
        size_t avail = m_BufferSize - m_BufferPos;
        _ASSERT(avail != 0);
        if (avail > nwrite) {
            avail = nwrite;
        }
        size_t advance = avail;
        if (src  &&  src != (const char*)(-1L)) {
            memcpy(m_Buffer + m_BufferPos, src, avail);
            size_t pad = ALIGN_SIZE(avail) - avail;
            memset(m_Buffer + m_BufferPos + avail, 0, pad);
            advance += pad;
            src += avail;
        }
        m_BufferPos += advance;
        if (m_BufferPos == m_BufferSize) {
            size_t nwritten = 0;
            do {
                streamsize xwritten = m_Stream->rdbuf()
                    ->sputn(             m_Buffer     + nwritten,
                            (streamsize)(m_BufferSize - nwritten));
                if (xwritten <= 0) {
                    int x_errno = errno;
                    m_Bad = true;
                    if (src != (const char*)(-1L)) {
                        TAR_THROW(this, eWrite,
                                  "Archive write failed" +s_OSReason(x_errno));
                    }
                    TAR_POST(84, Error,
                             "Archive write failed" + s_OSReason(x_errno));
                    return;
                }
                nwritten += xwritten;
            } while (nwritten < m_BufferSize);
            m_BufferPos = 0;
        }
        m_StreamPos += advance;
        nwrite -= avail;
    } while (nwrite);
}


// PAX (Portable Archive Interchange) extraction support

// Define bitmasks for extended numeric information
typedef enum {
    fPAXNone  = 0,
    fPAXMtime = 1 << 0,
    fPAXAtime = 1 << 1,
    fPAXCtime = 1 << 2,
    fPAXSize  = 1 << 3,
    fPAXUid   = 1 << 4,
    fPAXGid   = 1 << 5
} EPAXBit;
typedef unsigned int TPAXBits;  // Bitwise-OR of EPAXBit(s)


static bool s_ParsePAXInt(Uint8* valp, const char* str, size_t len, bool dot)
{
    if (!isdigit((unsigned char)(*str))) {
        return false;
    }
    char* e;
    errno = 0;
    unsigned long val = strtoul(str, &e, 10);
    if (errno  ||  (size_t)(e - str) > len  ||  (*e != '.'  &&  *e != '\n')) {
        return false;
    }
    if (*e++ == '.') {
        if (!dot) {
            return false;
        }
        if (!isdigit((unsigned char)(*e))) {
            return false;
        }
        _ASSERT(!errno);
        (void) strtoul(e, &e, 10);
        if (errno  ||  (size_t)(e - str) != len  ||  *e != '\n') {
            return false;
        }
    }
    *valp = val;
    return true;
}


static bool s_AllLowerCase(const char* str, size_t len)
{
    for (size_t i = 0;  i < len;  i++) {
        if (!islower((unsigned char) str[i]))
            return false;
    }
    return true;
}


CTar::EStatus CTar::x_ParsePAXHeader(const string& buffer)
{
    Uint8 mtime = 0, atime = 0, ctime = 0, size = 0, uid = 0, gid = 0;
    string path, linkpath, uname, gname;
    string* nodot = (string*)(-1);
    const struct SPAXParseTable {
        const char* key;
        Uint8*      val;  // non-null for numeric, else do as string
        string*     str;  // null for check only (numeric: non-null for no '.')
        EPAXBit     bit;  // for numerics only (NB: ePAXNone if check only)
    } parser[] = {
        { "mtime",    &mtime, 0,         fPAXMtime },  // numeric w/dot: assign
        { "atime",    &atime, 0,         fPAXAtime },
        { "ctime",    &ctime, 0,         fPAXCtime },
        { "size",     &size,  nodot,     fPAXSize  },  // num.-no-dot: assign
        { "uid",      &uid,   nodot,     fPAXUid   },
        { "gid",      &gid,   nodot,     fPAXGid   },
      /*{ "dummy",    &dummy, nodot,     fPAXNone  },*/// num.-no-dot: ck.only
        { "path",     0,      &path,     fPAXNone  },  // string: assign
        { "linkpath", 0,      &linkpath, fPAXNone  },
        { "uname",    0,      &uname,    fPAXNone  },
        { "gname",    0,      &gname,    fPAXNone  },
        { "comment",  0,      0,         fPAXNone  },  // string: check only
        { "charset",  0,      0,         fPAXNone  }   // string: check only
    };
    const char* str = buffer.c_str();
    TPAXBits parsed = fPAXNone;

    do {
        unsigned long len;
        size_t klen, vlen;
        char *k, *e, *v;

        errno = 0;
        if (!isdigit((unsigned char)(*str)) || !(e = (char*) strchr(str, '\n'))
            ||  !(len = strtoul(str, &k, 10))  ||  str + len - 1 != e
            ||  *k++ != ' '  ||  !(v = (char*) memchr(k, '=', (size_t)(e - k)))
            ||  !(klen = (size_t)(v++ - k))  ||  memchr(k, ' ', klen)
            ||  !(vlen = (size_t)(e - v))) {
            TAR_POST(74, Error, "Skipping malformed PAX header");
            return eFailure;
        }
        bool done = false;
        for (size_t n = 0;  n < sizeof(parser) / sizeof(parser[0]);  n++) {
            if (strlen(parser[n].key) == klen
                &&  strncmp(parser[n].key, k, klen) == 0) {
                if (!parser[n].val) {
                    if (parser[n].str)
                        parser[n].str->assign(v, vlen);
                } else if (!isdigit((unsigned char) v[0])
                           ||  !s_ParsePAXInt(parser[n].val, v, vlen,
                                              !parser[n].str ? true : false)) {
                    TAR_POST(75, Warning,
                             "Ignoring bad numeric '" + string(v, vlen)
                             + "' in PAX value '" + string(k, klen) + '\'');
                } else {
                    parsed |= parser[n].bit;
                }
                done = true;
                break;
            }
        }
        if (!done  &&  s_AllLowerCase(k, klen)/*&&  !memchr(k, '.', klen)*/) {
            TAR_POST(76, Warning,
                     "Ignoring unrecognized PAX value '"
                     + string(k, klen) + '\'');
        }
        str = ++e;
    } while (*str);

    m_Current.m_Name.swap(path);
    m_Current.m_LinkName.swap(linkpath);
    m_Current.m_UserName.swap(uname);
    m_Current.m_GroupName.swap(gname);
    m_Current.m_Stat.st_mtime = (time_t) mtime;
    m_Current.m_Stat.st_atime = (time_t) atime;
    m_Current.m_Stat.st_ctime = (time_t) ctime;
    m_Current.m_Stat.st_size  = (off_t)  size;
    m_Current.m_Stat.st_uid   = (uid_t)  uid;
    m_Current.m_Stat.st_gid   = (gid_t)  gid;
    m_Current.m_Pos = parsed;

    return eContinue;
}


static void s_Dump(const string& file, Uint8 pos, size_t recsize,
                   const string& entryname, const SHeader* h, ETar_Format fmt,
                   Uint8 datasize)
{
    EDiagSev level = SetDiagPostLevel(eDiag_Info);
    Uint8 blocks = BLOCK_OF(ALIGN_SIZE(datasize));
    ERR_POST(Info << '\n' + s_PositionAsString(file, pos, recsize, entryname)
             + s_DumpHeader(h, fmt) + '\n'
             + (blocks
                ? "Blocks of data: " + NStr::UInt8ToString(blocks) + '\n'
                : kEmptyStr));
    SetDiagPostLevel(level);
}


CTar::EStatus CTar::x_ReadEntryInfo(bool dump, bool pax)
{
    // Read block
    const TBlock* block;
    size_t nread = sizeof(block->buffer);
    _ASSERT(sizeof(*block) == BLOCK_SIZE/*== sizeof(block->buffer)*/);
    if (!(block = (const TBlock*) x_ReadArchive(nread))) {
        return eEOF;
    }
    if (nread != BLOCK_SIZE) {
        TAR_THROW(this, eRead,
                  "Unexpected EOF in archive");
    }
    const SHeader* h = &block->header;

    // Check header format
    ETar_Format fmt = eTar_Unknown;
    if (memcmp(h->magic, "ustar", 6) == 0) {
        fmt = pax ? eTar_Posix : eTar_Ustar;
    } else if (memcmp(h->magic, "ustar  ", 8) == 0) {
        // NB: Here, the magic protruded into the adjacent version field
        fmt = eTar_OldGNU;
    } else if (memcmp(h->magic, "\0\0\0\0\0", 6) == 0) {
        fmt = eTar_Legacy;
    } else {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Unrecognized format", h, fmt);
    }

    Uint8 val;
    // Get checksum from header
    if (!s_OctalToNum(val, h->checksum, sizeof(h->checksum))) {
        // We must allow all zero bytes here in case of pad/zero blocks
        for (size_t i = 0;  i < sizeof(block->buffer);  i++) {
            if (block->buffer[i]) {
                TAR_THROW_EX(this, eUnsupportedTarFormat,
                             "Bad checksum", h, fmt);
            }
        }
        m_StreamPos += BLOCK_SIZE;  // NB: nread
        return eZeroBlock;
    }
    int checksum = int(val);

    // Compute both signed and unsigned checksums (for compatibility)
    int ssum = 0;
    unsigned int usum = 0;
    const char* p = block->buffer;
    for (size_t i = 0;  i < sizeof(block->buffer);  i++)  {
        ssum +=                 *p;
        usum += (unsigned char)(*p);
        p++;
    }
    p = h->checksum;
    for (size_t j = 0;  j < sizeof(h->checksum);  j++) {
        ssum -=                 *p  - ' ';
        usum -= (unsigned char)(*p) - ' ';
        p++;
    }

    // Compare checksum(s)
    if (checksum != ssum   &&  (unsigned int) checksum != usum) {
        string message = "Header checksum failed";
        if (m_Flags & fDumpBlockHeaders) {
            message += ", expected ";
            if (usum != (unsigned int) ssum) {
                message += "either ";
            }
            if (usum > 7) {
                message += "0";
            }
            message += NStr::UIntToString(usum, 0, 8);
            if (usum != (unsigned int) ssum) {
                message += " or ";
                if ((unsigned int) ssum > 7) {
                    message += "0";
                }
                message += NStr::UIntToString((unsigned int) ssum, 0, 8);
            }
        }
        TAR_THROW_EX(this, eChecksum,
                     message, h, fmt);
    }

    // Set all info members now (thus, validating the header block)

    m_Current.m_HeaderSize = BLOCK_SIZE;

    // Name
    if (m_Current.GetName().empty()) {
        if ((fmt & eTar_Ustar)  &&  h->prefix[0]
            &&  tolower((unsigned char) h->typeflag[0]) != 'x') {
            m_Current.m_Name =
                CDirEntry::ConcatPath(string(h->prefix,
                                             s_Length(h->prefix,
                                                      sizeof(h->prefix))),
                                      string(h->name,
                                             s_Length(h->name,
                                                      sizeof(h->name))));
        } else {
            m_Current.m_Name.assign(h->name,
                                    s_Length(h->name, sizeof(h->name)));
        }
    }

    // Mode
    if (!s_OctalToNum(val, h->mode, sizeof(h->mode))) {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Bad entry mode", h, fmt);
    }
    m_Current.m_Stat.st_mode = (mode_t) val;

    // User Id
    if (!s_DecodeUint8(val, h->uid, sizeof(h->uid))) {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Bad user ID", h, fmt);
    }
    m_Current.m_Stat.st_uid = (uid_t) val;

    // Group Id
    if (!s_DecodeUint8(val, h->gid, sizeof(h->gid))) {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Bad group ID", h, fmt);
    }
    m_Current.m_Stat.st_gid = (gid_t) val;

    // Size
    if (!s_DecodeUint8(val, h->size, sizeof(h->size))) {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Bad entry size", h, fmt);
    }
    m_Current.m_Stat.st_size = (off_t) val;
    if (m_Current.GetSize() != val) {
        ERR_POST_ONCE(Critical << "CAUTION:"
                      " ***"
                      " This run-time may not support large TAR entries"
                      " (have you built it --with-lfs?)"
                      " ***");
    }

    // Modification time
    if (!s_OctalToNum(val, h->mtime, sizeof(h->mtime))) {
        TAR_THROW_EX(this, eUnsupportedTarFormat,
                     "Bad modification time", h, fmt);
    }
    m_Current.m_Stat.st_mtime = (time_t) val;

    if (fmt == eTar_OldGNU  ||  (fmt & eTar_Ustar)) {
        // User name
        m_Current.m_UserName.assign(h->uname,
                                    s_Length(h->uname, sizeof(h->uname)));
        // Group name
        m_Current.m_GroupName.assign(h->gname,
                                     s_Length(h->gname,sizeof(h->gname)));
    }

    if (fmt == eTar_OldGNU) {
        // Name prefix cannot be used because there are times, and other stuff
        // NB: times are valid for incremental archive only, so checks relaxed
        if (!s_OctalToNum(val, h->gt.atime, sizeof(h->gt.atime))) {
            if (memcchr(h->gt.atime, '\0', sizeof(h->gt.atime))) {
                TAR_THROW_EX(this, eUnsupportedTarFormat,
                             "Bad last access time", h, fmt);
            }
        } else {
            m_Current.m_Stat.st_atime = (time_t) val;
        }
        if (!s_OctalToNum(val, h->gt.ctime, sizeof(h->gt.ctime))) {
            if (memcchr(h->gt.ctime, '\0', sizeof(h->gt.ctime))) {
                TAR_THROW_EX(this, eUnsupportedTarFormat,
                             "Bad creation time", h, fmt);
            }
        } else {
            m_Current.m_Stat.st_ctime = (time_t) val;
        }
    }

    // Entry type
    switch (h->typeflag[0]) {
    case '7':
        ERR_POST_ONCE(Critical << "CAUTION:"
                      " *** Contiguous TAR entries processed as regular ones"
                      " ***");
        /*FALLTHRU*/
    case '0':
    case '\0':
        if (!(fmt & eTar_Ustar)  &&  fmt != eTar_OldGNU) {
            size_t namelen = s_Length(h->name, sizeof(h->name));
            if (namelen  &&  h->name[namelen - 1] == '/') {
                m_Current.m_Type = CTarEntryInfo::eDir;
                m_Current.m_Stat.st_size = 0;
                break;
            }
        }
        m_Current.m_Type = CTarEntryInfo::eFile;
        break;
    case '1':
    case '2':
        m_Current.m_Type = (h->typeflag[0] == '1'
                            ? CTarEntryInfo::eHardLink
                            : CTarEntryInfo::eSymLink);
        m_Current.m_LinkName.assign(h->linkname,
                                    s_Length(h->linkname,sizeof(h->linkname)));
        if (!m_Current.GetSize()) {
            break;
        }
        if (h->typeflag[0] != '1') {
            m_Current.m_Stat.st_size = 0;
        } else if (fmt != eTar_Posix) {
            TAR_POST(77, Warning,
                     "Non-zero hard-link size ("
                     + NStr::UInt8ToString(m_Current.GetSize())
                     + ") is ignored (non-PAX)");
            m_Current.m_Stat.st_size = 0;
        }
        break;
    case '3':
    case '4':
        m_Current.m_Type = (h->typeflag[0] == '3'
                            ? CTarEntryInfo::eCharDev
                            : CTarEntryInfo::eBlockDev);
        if (!s_OctalToNum(val, h->devminor, sizeof(h->devminor))) {
            TAR_THROW_EX(this, eUnsupportedTarFormat,
                         "Bad device minor number", h, fmt);
        }
        usum = (unsigned int) val; // set aside
        if (!s_OctalToNum(val, h->devmajor, sizeof(h->devmajor))) {
            TAR_THROW_EX(this, eUnsupportedTarFormat,
                         "Bad device major number", h, fmt);            
        }
#ifdef makedev
        m_Current.m_Stat.st_rdev = makedev((unsigned int) val, usum);
#else
        if (sizeof(int) >= 4  &&  sizeof(m_Current.m_Stat.st_rdev) >= 4) {
            *((unsigned int*) &m_Current.m_Stat.st_rdev) =
                (unsigned int)((val << 16) | usum);
        }
#endif // makedev
        m_Current.m_Stat.st_size = 0;
        break;
    case '5':
        m_Current.m_Type = CTarEntryInfo::eDir;
        m_Current.m_Stat.st_size = 0;
        break;
    case '6':
        m_Current.m_Type = CTarEntryInfo::ePipe;
        m_Current.m_Stat.st_size = 0;
        break;
    case 'x':
    case 'X':
    case 'K':
    case 'L':
        if ((tolower((unsigned char) h->typeflag[0]) == 'x'
             &&  (fmt & eTar_Ustar))  ||
            (tolower((unsigned char) h->typeflag[0]) != 'x'
             &&  fmt == eTar_OldGNU)) {
            // Assign actual type
            switch (h->typeflag[0]) {
            case 'x':
            case 'X':
                if (pax) {
                    TAR_POST(78, Warning,
                             "Double PAX header encountered,"
                             " archive may be corrupt");
                }
                fmt = eTar_Posix;  // upgrade
                m_Current.m_Type = CTarEntryInfo::ePAXHeader;
                break;
            case 'K':
                m_Current.m_Type = CTarEntryInfo::eGNULongLink;
                break;
            case 'L':
                m_Current.m_Type = CTarEntryInfo::eGNULongName;
                break;
            }
            // Dump header
            size_t hsize = (size_t) m_Current.GetSize();
            if (dump) {
                s_Dump(m_FileName, m_StreamPos, m_BufferSize,
                       m_Current.GetName(), h, fmt, hsize);
            }
            m_StreamPos += BLOCK_SIZE;  // NB: nread
            // Read in the extended information
            string buffer;
            while (hsize) {
                nread = hsize;
                const char* xbuf = x_ReadArchive(nread);
                if (!xbuf) {
                    TAR_THROW(this, eRead,
                              string("Unexpected EOF in ") +
                              (m_Current.GetType()
                               == CTarEntryInfo::ePAXHeader
                               ? "PAX header" :
                               m_Current.GetType()
                               == CTarEntryInfo::eGNULongName
                               ? "long name"
                               : "long link"));
                }
                buffer.append(xbuf, nread);
                m_StreamPos += ALIGN_SIZE(nread);
                hsize -= nread;
            }
            // Make sure there's no embedded '\0'(s)
            buffer.resize(strlen(buffer.c_str()));
            if (dump) {
                string what(m_Current.GetType() == CTarEntryInfo::ePAXHeader
                            ? "PAX header:\n" :
                            m_Current.GetType() == CTarEntryInfo::eGNULongName
                            ? "Long name:      \""
                            : "Long link name: \"");
                EDiagSev level = SetDiagPostLevel(eDiag_Info);
                ERR_POST(Info << Message << '\n' + what
                         + NStr::PrintableString(buffer,
                                                 m_Current.GetType()
                                                 == CTarEntryInfo::ePAXHeader
                                                 ? NStr::fNewLine_Passthru
                                                 : NStr::fNewLine_Quote) +
                         (m_Current.GetType() == CTarEntryInfo::ePAXHeader ?
                          buffer.size()  &&  buffer[buffer.size()-1] == '\n'
                          ? kEmptyStr : "\n"
                          : "\"\n"));
                SetDiagPostLevel(level);
            }
            // Reset size because the data blocks have been all read
            hsize = (size_t) m_Current.GetSize();
            m_Current.m_HeaderSize += ALIGN_SIZE(hsize);
            m_Current.m_Stat.st_size = 0;
            if (!hsize  ||  !buffer.size()) {
                TAR_POST(79, Error,
                         "Skipping " + string(hsize ? "empty" : "zero-sized")
                         + " extended header");
                return eFailure;
            }
            if (m_Current.GetType() == CTarEntryInfo::ePAXHeader) {
                if (hsize != buffer.size()) {
                    TAR_POST(80, Error,
                             "Skipping truncated ("
                             + NStr::UInt8ToString((Uint8) hsize) + "->"
                             + NStr::UInt8ToString((Uint8) buffer.size())
                             + ") PAX header");
                    return eFailure;
                }
                return x_ParsePAXHeader(buffer);
            }
            if (m_Current.GetType() == CTarEntryInfo::eGNULongName) {
                m_Current.m_Name.swap(buffer);
            } else {
                m_Current.m_LinkName.swap(buffer);
            }
            return eContinue;
        }
        /*FALLTHRU*/
    default:
        m_Current.m_Type = CTarEntryInfo::eUnknown;
        break;
    }

    if (dump) {
        s_Dump(m_FileName, m_StreamPos, m_BufferSize,
               m_Current.GetName(), h, fmt, m_Current.GetSize());
    }
    m_StreamPos += BLOCK_SIZE;  // NB: nread

    return eSuccess;
}


void CTar::x_WriteEntryInfo(const string& name)
{
    // Prepare block info
    TBlock block;
    _ASSERT(sizeof(block) == BLOCK_SIZE/*== sizeof(block.buffer)*/);
    memset(block.buffer, 0, sizeof(block.buffer));
    SHeader* h = &block.header;

    CTarEntryInfo::EType type = m_Current.GetType();

    // Name(s) ('\0'-terminated if fit entirely, otherwise not)
    if (!x_PackName(h, m_Current, false)) {
        TAR_THROW(this, eNameTooLong,
                  "Name '" + m_Current.GetName() + "' too long in"
                  " entry '" + name + '\'');
    }
    if (type == CTarEntryInfo::eSymLink  &&  !x_PackName(h, m_Current, true)) {
        TAR_THROW(this, eNameTooLong,
                  "Link '" + m_Current.GetLinkName() + "' too long in"
                  " entry '" + name + '\'');
    }

    /* NOTE:  Although some sources on the Internet indicate that
     * all but size, mtime and version numeric fields are '\0'-terminated,
     * we could not confirm that with existing tar programs, all of
     * which we saw using either '\0' or ' '-terminated values in both size
     * and mtime fields.  For ustar archive we have found a document
     * that definitively tells that _all_ numeric fields are '\0'-terminated,
     * and can keep maximum sizeof(field)-1 octal digits.  We follow it here.
     * However, GNU and ustar checksums seem to be different indeed,
     * so we don't use a trailing space for ustar, but for GNU only.
     */

    // Mode
    if (!s_NumToOctal(m_Current.GetMode(), h->mode, sizeof(h->mode) - 1)) {
        TAR_THROW(this, eMemory,
                  "Cannot store file mode");
    }

    // Update format as we go
    ETar_Format fmt = eTar_Ustar;
    int ok;

    // User ID
    ok = s_EncodeUint8(m_Current.GetUserId(), h->uid, sizeof(h->uid) - 1);
    if (!ok) {
        TAR_THROW(this, eMemory,
                  "Cannot store user ID");
    }
    if (ok < 0) {
        fmt = eTar_OldGNU;
    }

    // Group ID
    ok = s_EncodeUint8(m_Current.GetGroupId(), h->gid, sizeof(h->gid) - 1);
    if (!ok) {
        TAR_THROW(this, eMemory,
                  "Cannot store group ID");
    }
    if (ok < 0) {
        fmt = eTar_OldGNU;
    }

    // Size
    _ASSERT(type == CTarEntryInfo::eFile  ||  m_Current.GetSize() == 0);
    ok = s_EncodeUint8(m_Current.GetSize(), h->size, sizeof(h->size) - 1);
    if (!ok) {
        TAR_THROW(this, eMemory,
                  "Cannot store file size");
    }
    if (ok < 0) {
        fmt = eTar_OldGNU;
    }

    if (fmt != eTar_Ustar  &&  h->prefix[0]) {
        // cannot downgrade to reflect encoding
        fmt  = eTar_Ustar;
    }

    // Modification time
    if (!s_NumToOctal(m_Current.GetModificationTime(),
                      h->mtime, sizeof(h->mtime) - 1)) {
        TAR_THROW(this, eMemory,
                  "Cannot store modification time");
    }

    bool device = false;
    // Type (GNU extension for SymLink)
    switch (type) {
    case CTarEntryInfo::eFile:
        h->typeflag[0] = '0';
        break;
    case CTarEntryInfo::eSymLink:
        h->typeflag[0] = '2';
        break;
    case CTarEntryInfo::eCharDev:
    case CTarEntryInfo::eBlockDev:
        h->typeflag[0] = type == CTarEntryInfo::eCharDev ? '3' : '4';
        if (!s_NumToOctal(m_Current.GetMajor(),
                          h->devmajor, sizeof(h->devmajor) - 1)) {
            TAR_THROW(this, eMemory,
                      "Cannot store major number");
        }
        if (!s_NumToOctal(m_Current.GetMinor(),
                          h->devminor, sizeof(h->devminor) - 1)) {
            TAR_THROW(this, eMemory,
                      "Cannot store minor number");
        }
        device = true;
        break;
    case CTarEntryInfo::eDir:
        h->typeflag[0] = '5';
        break;
    case CTarEntryInfo::ePipe:
        h->typeflag[0] = '6';
        break;
    default:
        TAR_THROW(this, eUnsupportedEntryType,
                  "Don't know how to store entry '" + name + "' w/type #" +
                  NStr::IntToString(int(type)) + " into archive: "
                  "Internal error, please report!");
    }

    // User and group
    const string& usr = m_Current.GetUserName();
    size_t len = usr.length();
    if (len < sizeof(h->uname)) {
        memcpy(h->uname, usr.c_str(), len);
    }
    const string& grp = m_Current.GetGroupName();
    len = grp.length();
    if (len < sizeof(h->gname)) {
        memcpy(h->gname, grp.c_str(), len);
    }

    // Device nos to complete the ustar header protocol (all fields ok)
    if (!device  &&  fmt != eTar_OldGNU) {
        s_NumToOctal(0, h->devmajor, sizeof(h->devmajor) - 1);
        s_NumToOctal(0, h->devminor, sizeof(h->devminor) - 1);
    }

    if (fmt != eTar_OldGNU) {
        // Magic
        strcpy(h->magic,   "ustar");
        // Version (EXCEPTION:  not '\0' terminated)
        memcpy(h->version, "00", 2);
    } else {
        // NB: Old GNU magic protrudes into adjacent version field
        memcpy(h->magic,   "ustar  ", 8); // 2 spaces and '\0'-terminated
    }

    // NCBI signature if allowed
    if (!(m_Flags & fStandardHeaderOnly)) {
        _ASSERT(sizeof(block.header) + 4 < sizeof(block.buffer));
        memcpy(block.buffer + sizeof(block) - 4, "NCBI", 4);
    }

    // Final step: checksumming
    if (!s_TarChecksum(&block, fmt == eTar_OldGNU ? true : false)) {
        TAR_THROW(this, eMemory,
                  "Cannot store checksum");
    }

    // Write header
    x_WriteArchive(sizeof(block.buffer), block.buffer);
    m_Current.m_HeaderSize = (streamsize)(m_StreamPos - m_Current.m_Pos);

    Checkpoint(m_Current, true/*write*/);
}


bool CTar::x_PackName(SHeader* h, const CTarEntryInfo& info, bool link)
{
    char*      storage = link ? h->linkname         : h->name;
    size_t        size = link ? sizeof(h->linkname) : sizeof(h->name);
    const string& name = link ? info.GetLinkName()  : info.GetName();
    size_t         len = name.length();
    const char*    src = name.c_str();

    if (len <= size) {
        // Name fits!
        memcpy(storage, src, len);
        return true;
    }

    if (!link  &&  len <= sizeof(h->prefix) + 1 + sizeof(h->name)) {
        // Try to split the long name into a prefix and a short name (POSIX)
        size_t i = len;
        if (i > sizeof(h->prefix)) {
            i = sizeof(h->prefix);
        }
        while (i > 0  &&  src[--i] != '/');
        if (i  &&  len - i <= sizeof(h->name) + 1) {
            memcpy(h->prefix, src,         i);
            memcpy(h->name,   src + i + 1, len - i - 1);
            return true;
        }
    }

    // Still, store the initial part in the original header
    memcpy(storage, src, size);

    // Prepare extended block header with the long name info (old GNU style)
    _ASSERT(!OFFSET_OF(m_BufferPos)  &&  m_BufferPos < m_BufferSize);
    TBlock* block = (TBlock*)(m_Buffer + m_BufferPos);
    memset(block->buffer, 0, sizeof(block->buffer));
    h = &block->header;

    // See above for comments about header filling
    len++; // write terminating '\0' as it can always be made to fit in
    strcpy(h->name, "././@LongLink");
    s_NumToOctal(0,         h->mode,  sizeof(h->mode) - 1);
    s_NumToOctal(0,         h->uid,   sizeof(h->uid)  - 1);
    s_NumToOctal(0,         h->gid,   sizeof(h->gid)  - 1);
    if (!s_EncodeUint8(len, h->size,  sizeof(h->size) - 1)) {
        return false;
    }
    s_NumToOctal(0,         h->mtime, sizeof(h->mtime)- 1);
    h->typeflag[0] = link ? 'K' : 'L';

    // NB: Old GNU magic protrudes into adjacent version field
    memcpy(h->magic, "ustar  ", 8); // 2 spaces and '\0'-terminated

    s_TarChecksum(block, true);

    // Write the header
    x_WriteArchive(sizeof(block->buffer));

    // Store the full name in the extended block
    AutoPtr< char, ArrayDeleter<char> > buf_ptr(new char[len]);
    storage = buf_ptr.get();
    memcpy(storage, src, len);

    // Write the extended block (will be aligned as necessary)
    x_WriteArchive(len, storage);

    return true;
}


void CTar::x_Backspace(EAction action, Uint8 blocks)
{
    _ASSERT(SIZE_OF(blocks) <= m_StreamPos);
    m_Current.m_Name.erase();
    if (!blocks  ||  (action != eAppend  &&  action != eUpdate)) {
        return;
    }
    if (!m_FileStream) {
        TAR_POST(4, Warning,
                 "In-stream update may result in gapped tar archive");
        return;
    }

    // Assertion:  there's data (perhaps 0's) in the file, backup is possible.
    _ASSERT(!OFFSET_OF(m_StreamPos));
    Uint8 gap = SIZE_OF(blocks);
    m_StreamPos -= gap;
    if (m_BufferPos == 0) {
        m_BufferPos += m_BufferSize;
    }
    CT_POS_TYPE rec  = (CT_OFF_TYPE)(m_StreamPos / m_BufferSize);
    size_t      off  = (size_t)     (m_StreamPos % m_BufferSize);
    if (gap > m_BufferPos) {
        size_t temp  = 0;
        m_BufferPos  = 0;
        // Refetch the block
        if (!m_FileStream->seekg(rec * m_BufferSize)  ||
            !x_ReadArchive(temp = BLOCK_SIZE)) {
            TAR_POST(85, Error,
                     "Archive backspace error in record reget"
                     + string(temp ? " (EOF)" : ""));
        }
        m_BufferPos  = off;
    } else {
        m_BufferPos -= (size_t) gap;
    }

    // Always set put position here
    if (!m_FileStream->seekp(rec * m_BufferSize)) {
        TAR_POST(86, Error,
                 "Archive backspace error in record reset");
    }
}


auto_ptr<CTar::TEntries> CTar::x_ReadAndProcess(EAction action)
{
    auto_ptr<TEntries> done(new TEntries);
    size_t zeroblock_count = 0;
    Uint8 pos = m_StreamPos;
    CTarEntryInfo xinfo;

    for (;;) {
        // Next block is supposed to be a header
        m_Current = CTarEntryInfo(pos);
        m_Current.m_Name = xinfo.GetName();
        EStatus status = x_ReadEntryInfo
            (action == eTest  &&  (m_Flags & fDumpBlockHeaders),
             xinfo.GetType() == CTarEntryInfo::ePAXHeader);
        switch (status) {
        case eFailure:
        case eSuccess:
        case eContinue:
            if (zeroblock_count  &&  !(m_Flags & fIgnoreZeroBlocks)) {
                TAR_POST(5, Error, "Interspersing single zero block ignored");
            }
            break;

        case eZeroBlock:
            zeroblock_count++;
            if ((m_Flags & fIgnoreZeroBlocks)  ||  zeroblock_count < 2) {
                if (xinfo.GetType() == CTarEntryInfo::eUnknown) {
                    // Not yet reading an entry -- advance
                    pos += BLOCK_SIZE;
                }
                continue;
            }
            // Two zero blocks -> eEOF
            /*FALLTHRU*/

        case eEOF:
            if (xinfo.GetType() != CTarEntryInfo::eUnknown) {
                TAR_POST(6, Error, "Orphaned extended information ignored");
            }
            x_Backspace(action, zeroblock_count);
            return done;
        }
        zeroblock_count = 0;

        //
        // Process entry
        //
        if (status == eContinue) {
            // Extended header information has just been read in
            xinfo.m_HeaderSize += m_Current.m_HeaderSize;

            switch (m_Current.GetType()) {
            case CTarEntryInfo::ePAXHeader:
                if (xinfo.GetType() != CTarEntryInfo::eUnknown) {
                    TAR_POST(7, Error, "Unused extended header replaced");
                }
                xinfo.m_Type = CTarEntryInfo::ePAXHeader;
                xinfo.m_Name.swap(m_Current.m_Name);
                xinfo.m_LinkName.swap(m_Current.m_LinkName);
                xinfo.m_UserName.swap(m_Current.m_UserName);
                xinfo.m_GroupName.swap(m_Current.m_GroupName);
                xinfo.m_Stat = m_Current.m_Stat;
                xinfo.m_Pos  = m_Current.m_Pos;  // NB: Parsed mask, not pos!
                break;

            case CTarEntryInfo::eGNULongName:
                if (xinfo.GetType() == CTarEntryInfo::ePAXHeader
                    ||  !xinfo.GetName().empty()) {
                    TAR_POST(8, Error,
                             "Unused long name '" + xinfo.GetName()
                             + "' replaced");
                }
                // Latch next long name here then just skip
                xinfo.m_Type = CTarEntryInfo::eGNULongName;
                xinfo.m_Name.swap(m_Current.m_Name);
                break;

            case CTarEntryInfo::eGNULongLink:
                if (xinfo.GetType() == CTarEntryInfo::ePAXHeader
                    ||  !xinfo.GetLinkName().empty()) {
                    TAR_POST(9, Error,
                             "Unused long link '" + xinfo.GetLinkName()
                             + "' replaced");
                }
                // Latch next long link here then just skip
                xinfo.m_Type = CTarEntryInfo::eGNULongLink;
                xinfo.m_LinkName.swap(m_Current.m_LinkName);
                break;

            default:
                NCBI_THROW(CCoreException, eCore, "Internal error");
            }
            continue;
        }

        // Fixup current 'info' with extended information obtained previously
        m_Current.m_HeaderSize += xinfo.m_HeaderSize;
        xinfo.m_HeaderSize = 0;
        if (!xinfo.GetName().empty()) {
            xinfo.m_Name.swap(m_Current.m_Name);
            xinfo.m_Name.erase();
        }
        if (!xinfo.GetLinkName().empty()) {
            xinfo.m_LinkName.swap(m_Current.m_LinkName);
            xinfo.m_LinkName.erase();
        }
        if (xinfo.GetType() == CTarEntryInfo::ePAXHeader) {
            TPAXBits parsed = (TPAXBits) xinfo.m_Pos;
            if (!xinfo.GetUserName().empty()) {
                xinfo.m_UserName.swap(m_Current.m_UserName);
                xinfo.m_UserName.erase();
            }
            if (!xinfo.GetGroupName().empty()) {
                xinfo.m_GroupName.swap(m_Current.m_GroupName);
                xinfo.m_GroupName.erase();
            }
            if (parsed & fPAXMtime) {
                m_Current.m_Stat.st_mtime = xinfo.m_Stat.st_mtime;
            }
            if (parsed & fPAXAtime) {
                m_Current.m_Stat.st_atime = xinfo.m_Stat.st_atime;
            }
            if (parsed & fPAXCtime) {
                m_Current.m_Stat.st_ctime = xinfo.m_Stat.st_ctime;
            }
            if (parsed & fPAXSize) {
                m_Current.m_Stat.st_size = xinfo.m_Stat.st_size;
            }
            if (parsed & fPAXUid) {
                m_Current.m_Stat.st_uid = xinfo.m_Stat.st_uid;
            }
            if (parsed & fPAXGid) {
                m_Current.m_Stat.st_gid = xinfo.m_Stat.st_gid;
            }
        }
        xinfo.m_Type = CTarEntryInfo::eUnknown;
        _ASSERT(status == eFailure  ||  status == eSuccess);

        if (!Checkpoint(m_Current, false/*read*/))
            status = eFailure;

        // Match file name with the set of masks
        bool match = (status == eFailure ? false :
                      m_Mask  &&  (action == eList     ||
                                   action == eExtract  ||
                                   action == eInternal)
                      ? m_Mask->Match(m_Current.GetName(),
                                      m_Flags & fMaskNocase
                                      ? NStr::eNocase
                                      : NStr::eCase)
                      : true);

        // NB: match is 'false' when processing a failing entry
        if ((match  &&  action == eInternal)
            ||  x_ProcessEntry(match  &&  action == eExtract, done.get())
            ||  (match  &&  !(int(action) & ((eExtract | eInternal) & ~eRW)))){
            _ASSERT(status == eSuccess);
            done->push_back(m_Current);
            if (action == eInternal) {
                break;
            }
        }

        pos = m_StreamPos;
    }

    return done;
}


struct CTmpDirEntryDeleter {
    static void Delete(CDirEntry* entry)
    { entry->Remove(); delete entry; }
};


bool CTar::x_ProcessEntry(bool extract, const CTar::TEntries* done)
{
    Uint8                size = m_Current.GetSize();
    CTarEntryInfo::EType type = m_Current.GetType();

    if (extract) {
        // Destination for extraction
        auto_ptr<CDirEntry> dst
            (CDirEntry::CreateObject(CDirEntry::EType(type),
                                     CDirEntry::NormalizePath
                                     (CDirEntry::ConcatPath
                                      (m_BaseDir, m_Current.GetName()))));
        // Source for extraction
        auto_ptr<CDirEntry> src;
        // Direntry pending removal
        AutoPtr<CDirEntry, CTmpDirEntryDeleter> pending;

        // Dereference sym.link if requested
        if (type != CTarEntryInfo::eSymLink  &&
            type != CTarEntryInfo::eHardLink  &&  (m_Flags & fFollowLinks)) {
            dst->DereferenceLink();
        }

        // Actual type in file system (if exists)
        CDirEntry::EType dst_type = dst->GetType();

        // Look if extraction is allowed (when the destination exists)
        if (dst_type != CDirEntry::eUnknown) {
            bool found = false;  // check if ours (prev. revision extracted)
            if (done) {
                ITERATE(TEntries, i, *done) {
                    if (i->GetName() == m_Current.GetName()  &&
                        i->GetType() == m_Current.GetType()) {
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                // Can overwrite it?
                if (!(m_Flags & fOverwrite)) {
                    // File already exists, and cannot be changed
                    extract = false;
                } else { // The fOverwrite flag is set
                    // Can update?
                    if ((m_Flags & fUpdate) == fUpdate
                        &&  type != CTarEntryInfo::eDir) {
                        // Update directories always, because the archive can
                        // contain other subtree of this existing directory.
                        time_t dst_time;
                        // Make sure that dst is not older than the entry
                        if (dst->GetTimeT(&dst_time)
                            &&  dst_time >= m_Current.GetModificationTime()) {
                            extract = false;
                        }
                    }
                    // Have equal types?
                    if (extract  &&  (m_Flags & fEqualTypes)) {
                        if (type == CTarEntryInfo::eHardLink) {
                            src.reset
                                (new CDirEntry(CDirEntry::NormalizePath
                                               (CDirEntry::ConcatPath
                                                (m_BaseDir,
                                                 m_Current.GetLinkName()))));
                            if (dst_type != src->GetType()) {
                                extract = false;
                            }
                        } else if (dst_type != CDirEntry::EType(type)) {
                            extract = false;
                        }
                    }
                }
            }
            if (extract) {
                if (!found  &&  (m_Flags & fBackup) == fBackup) {
                    // Need to backup the existing destination?
                    CDirEntry tmp(*dst);
                    if (!tmp.Backup(kEmptyStr, CDirEntry::eBackup_Rename)) {
                        TAR_THROW(this, eBackup,
                                  "Failed to backup '" + dst->GetPath() +'\'');
                    }
                } else if (type != CTarEntryInfo::eDir) {
                    // Do removal safely -- until extraction is confirmed
                    CDirEntry tmp(*dst);
                    pending.reset(new CDirEntry(CDirEntry::GetTmpNameEx
                                                (dst->GetDir(), "xNCBItArX")));
                    errno = 0;
                    if (!tmp.Rename(pending->GetPath())  ||  dst->Exists()) {
                        // Security concern:  do not attempt data extraction
                        // into special files etc., which can harm the system.
                        int x_errno = errno ? errno : EEXIST;
                        TAR_THROW(this, eWrite,
                                  "Cannot extract '" + dst->GetPath() + '\''
                                  + s_OSReason(x_errno));
                    }
                }
            }
        }
        if (extract) {
#ifdef NCBI_OS_UNIX
            mode_t u;
            u = umask(0);
            umask(u & 077);
            try {
#endif // NCBI_OS_UNIX
                extract = x_ExtractEntry(size, dst.get(), src.get());
#ifdef NCBI_OS_UNIX
            } catch (...) {
                umask(u);
                throw;
            }
            umask(u);
#endif // NCBI_OS_UNIX
            if (!extract  &&  pending.get()/*NB: not dir*/) {
                dst->Remove();
                // Undo delete
                CDirEntry tmp(*pending);
                if (!tmp.Rename(dst->GetPath())) {
                    int x_errno = errno;
                    TAR_THROW(this, eWrite,
                              "Cannot restore '" + dst->GetPath()
                              + "' back in place" + s_OSReason(x_errno));
                }
            }
        }
    }

    x_Skip(BLOCK_OF(ALIGN_SIZE(size)));

    return extract;
}


void CTar::x_Skip(Uint8 blocks)
{
    while (blocks) {
#ifndef NCBI_COMPILER_WORKSHOP
        // RogueWave RTL is buggy in seeking pipes -- it clobbers
        // (discards) streambuf data instead of leaving it alone..
        if (!(m_Flags & fSlowSkipWithRead)
            &&  !m_BufferPos  &&  blocks >= BLOCK_OF(m_BufferSize)) {
            CT_OFF_TYPE fskip =
                (CT_OFF_TYPE)(blocks / BLOCK_OF(m_BufferSize) * m_BufferSize);
            _ASSERT(ALIGN_SIZE(fskip) == fskip);
            if (m_Stream->rdbuf()->PUBSEEKOFF(fskip, IOS_BASE::cur)
                != (CT_POS_TYPE)((CT_OFF_TYPE)(-1))) {
                blocks      -= BLOCK_OF(fskip);
                m_StreamPos +=          fskip;
                continue;
            }
            if (m_FileStream) {
                TAR_POST(2, Warning,
                         "Cannot fast skip in file archive '" +
                         m_FileName + "', reverting to slow skip");
            }
            m_Flags |= fSlowSkipWithRead;
        }
#endif // NCBI_COMPILER_WORKSHOP
        size_t nskip = (blocks < BLOCK_OF(m_BufferSize)
                        ? (size_t) SIZE_OF(blocks)
                        : m_BufferSize);
        _ASSERT(ALIGN_SIZE(nskip) == nskip);
        if (!x_ReadArchive(nskip)) {
            TAR_THROW(this, eRead,
                      "Archive skip failed (EOF)");
        }
        nskip        = ALIGN_SIZE(nskip);
        blocks      -= BLOCK_OF(nskip);
        m_StreamPos +=          nskip;
    }
}


bool CTar::x_ExtractEntry(Uint8& size,
                          const CDirEntry* dst, const CDirEntry* src)
{
    CTarEntryInfo::EType type = m_Current.GetType();
    auto_ptr<CDirEntry> src_ptr;  // deleter
    bool result = true;  // assume best

    if (type == CTarEntryInfo::eUnknown  &&  !(m_Flags & fSkipUnsupported)) {
        // Conform to POSIX-mandated behavior to extract as files
        type = CTarEntryInfo::eFile;
    }
    switch (type) {
    case CTarEntryInfo::eHardLink:
    case CTarEntryInfo::eFile:
        {{
            // Create base directory
            CDir dir(dst->GetDir());
            if (!dir.CreatePath()) {
                int x_errno = errno;
                TAR_THROW(this, eCreate,
                          "Cannot create directory '" + dir.GetPath() + '\''
                          + s_OSReason(x_errno));
            }

            if (type != CTarEntryInfo::eFile) {
                if (!src) {
                    src_ptr.reset(new CDirEntry(CDirEntry::NormalizePath
                                                (CDirEntry::ConcatPath
                                                 (m_BaseDir,
                                                  m_Current.GetLinkName()))));
                    src = src_ptr.get();
                }
                if (src->GetType() == CDirEntry::eUnknown  &&  size) {
                    // Looks like a dangling hard link but luckily we have
                    // the actual file data (POSIX extension) so use it here.
                    type = CTarEntryInfo::eFile;
                }
            }

            if (type == CTarEntryInfo::eFile) {
                // Create the file
                // FIXME:  Switch to CFileIO eventually to bypass
                // ofstream obscurity w.r.t. errors, extra buffering etc.
                // FIXME:  Should the file name match an existing device (or
                // a terminal, in particular), things may go really ugly here.
                ofstream ofs(dst->GetPath().c_str(),
                             IOS_BASE::out    |
                             IOS_BASE::binary |
                             IOS_BASE::trunc);
                if (!ofs) {
                    int x_errno = errno;
                    TAR_THROW(this, eCreate,
                              "Cannot create file '" + dst->GetPath() + '\''
                              + s_OSReason(x_errno));
                }

                while (size) {
                    // Read from the archive
                    size_t nread = size < m_BufferSize
                        ? (size_t) size : m_BufferSize;
                    const char* xbuf = x_ReadArchive(nread);
                    if (!xbuf) {
                        TAR_THROW(this, eRead,
                                  "Unexpected EOF");
                    }
                    // Write file to disk
                    if (!ofs.write(xbuf, (streamsize) nread)) {
                        int x_errno = errno;
                        TAR_THROW(this, eWrite,
                                  "Error writing file '" +dst->GetPath()+ '\''
                                  + s_OSReason(x_errno));
                    }
                    m_StreamPos += ALIGN_SIZE(nread);
                    size -= nread;
                }

                _ASSERT(ofs.good());
                ofs.close();
            } else {
                _ASSERT(src);
#ifdef NCBI_OS_UNIX
                if (link(src->GetPath().c_str(),
                         dst->GetPath().c_str()) == 0) {
                    if (m_Flags & fPreserveAll) {
                        x_RestoreAttrs(m_Current, dst);
                    }
                    break;
                }
                int x_errno = errno;
                TAR_POST(10, Warning,
                         "Cannot hard-link '" + src->GetPath()
                         + "' and '" + dst->GetPath() + '\''
                         + s_OSReason(x_errno) + ", trying to copy");
#endif // NCBI_OS_UNIX
                if (!src->Copy(dst->GetPath(),
                               CDirEntry::fCF_Overwrite |
                               CDirEntry::fCF_PreserveAll)) {
                    TAR_POST(11, Error,
                             "Cannot hard-link '" + src->GetPath()
                             + "' and '" + dst->GetPath() + "\' via copy");
                    result = false;
                    break;
                }
            }

            // Restore attributes
            if (m_Flags & fPreserveAll) {
                x_RestoreAttrs(m_Current, dst);
            }
        }}
        break;

    case CTarEntryInfo::eDir:
        if (!CDir(dst->GetPath()).CreatePath()) {
            int x_errno = errno;
            TAR_THROW(this, eCreate,
                      "Cannot create directory '" + dst->GetPath() + '\''
                      + s_OSReason(x_errno));
        }
        // Attributes for a directory must be set only when all
        // its files have been already extracted.
        _ASSERT(size == 0);
        break;

    case CTarEntryInfo::eSymLink:
        {{
            CSymLink symlink(dst->GetPath());
            if (!symlink.Create(m_Current.GetLinkName())) {
                int x_errno = errno;
                TAR_POST(12, Error,
                         "Cannot create symlink '" + dst->GetPath()
                         + "' -> '" + m_Current.GetLinkName() + '\''
                         + s_OSReason(x_errno));
                result = false;
            }
            _ASSERT(size == 0);
        }}
        break;

    case CTarEntryInfo::ePAXHeader:
    case CTarEntryInfo::eGNULongName:
    case CTarEntryInfo::eGNULongLink:
        // Extended headers should have already been processed and not be here
        _TROUBLE;
        /*FALLTHRU*/

    case CTarEntryInfo::ePipe:
        {{
            _ASSERT(size == 0);
#ifdef NCBI_OS_UNIX
            mode_t u = umask(0);
            if (mkfifo(dst->GetPath().c_str(), m_Current.GetMode()) != 0) {
                result = false;
            }
            umask(u);  // NB: always succeeds and does not change errno
            if (result) {
                break;
            }
            string reason = s_OSReason(errno);
#else
            string reason = ": Feature not supported by host OS";
            result = false;
#endif // NCBI_OS_UNIX
            TAR_POST(81, Error,
                     "Cannot create FIFO '" + dst->GetPath() + '\'' + reason);
        }}
        break;

    case CTarEntryInfo::eCharDev:
    case CTarEntryInfo::eBlockDev:
        {{
            _ASSERT(size == 0);
#ifdef NCBI_OS_UNIX
            mode_t u = umask(0);
            mode_t m = (m_Current.GetMode() |
                        (type == CTarEntryInfo::eCharDev ? S_IFCHR : S_IFBLK));
            if (mknod(dst->GetPath().c_str(), m, m_Current.m_Stat.st_rdev)) {
                result = false;

            }
            umask(u);  // NB: always succeeds and does not clobber errno
            if (result) {
                break;
            }
            string reason = s_OSReason(errno);
#else
            string reason = ": Feature not supported by host OS";
            result = false;
#endif // NCBI_OS_UNIX
            TAR_POST(82, Error,
                     "Cannot create " +
                     string(type == CTarEntryInfo::eCharDev
                            ? "character" : "block")
                     + " special '" + dst->GetPath() + '\'' + reason);
        }}
        break;

    default:
        TAR_POST(13, Warning,
                 "Skipping unsupported entry '" + m_Current.GetName()
                 + "' w/type #" + NStr::IntToString(int(type)));
        result = false;
        break;
    }

    return result;
}


void CTar::x_RestoreAttrs(const CTarEntryInfo& info,
                          const CDirEntry*     dst) const
{
    auto_ptr<CDirEntry> dst_ptr;  // deleter
    if (!dst) {
        dst_ptr.reset(CDirEntry::CreateObject
                      (CDirEntry::EType(info.GetType()),
                       CDirEntry::NormalizePath
                       (CDirEntry::ConcatPath
                        (m_BaseDir, info.GetName()))));
        dst = dst_ptr.get();
    }

    // Date/time.
    // Set the time before permissions because on some platforms
    // this setting can also affect file permissions.
    if (m_Flags & fPreserveTime) {
        time_t modification(info.GetModificationTime());
        time_t last_access(info.GetLastAccessTime());
        time_t creation(info.GetCreationTime());
        if (!dst->SetTimeT(&modification, &last_access, &creation)) {
            int x_errno = errno;
            TAR_THROW(this, eRestoreAttrs,
                      "Cannot restore date/time for '" + dst->GetPath() + '\''
                      + s_OSReason(x_errno));
        }
    }

    // Owner.
    // This must precede changing permissions because on some
    // systems chown() clears the set[ug]id bits for non-superusers
    // thus resulting in incorrect permissions.
    if (m_Flags & fPreserveOwner) {
        // 2-tier trial:  first using the names, then using numeric IDs.
        // Note that it is often impossible to restore the original owner
        // without super-user rights so no error checking is done here.
        if (!dst->SetOwner(info.GetUserName(),
                           info.GetGroupName(),
                           eIgnoreLinks)) {
            dst->SetOwner(NStr::IntToString(info.GetUserId()),
                          NStr::IntToString(info.GetGroupId()),
                          eIgnoreLinks);
        }
    }

    // Mode.
    // Set them last.
    if ((m_Flags & fPreserveMode)
        &&  info.GetType() != CTarEntryInfo::ePipe
        &&  info.GetType() != CTarEntryInfo::eCharDev
        &&  info.GetType() != CTarEntryInfo::eBlockDev) {
        bool failed = false;
#ifdef NCBI_OS_UNIX
        // We cannot change permissions for sym.links because lchmod()
        // is not portable and is not implemented on majority of platforms.
        if (info.GetType() != CTarEntryInfo::eSymLink) {
            // Use raw mode here to restore most of the bits
            mode_t mode = s_TarToMode(info.m_Stat.st_mode);
            if (chmod(dst->GetPath().c_str(), mode) != 0) {
                // May fail due to setuid/setgid bits -- strip'em and try again
                if (mode &   (S_ISUID | S_ISGID)) {
                    mode &= ~(S_ISUID | S_ISGID);
                    failed = chmod(dst->GetPath().c_str(), mode) != 0;
                } else {
                    failed = true;
                }
            }
        }
#else
        CDirEntry::TMode user, group, other;
        CDirEntry::TSpecialModeBits special_bits;
        info.GetMode(&user, &group, &other, &special_bits);
        failed = !dst->SetMode(user, group, other, special_bits);
#endif // NCBI_OS_UNIX
        if (failed) {
            int x_errno = errno;
            TAR_THROW(this, eRestoreAttrs,
                      "Cannot restore mode bits for '" + dst->GetPath() + '\''
                      + s_OSReason(x_errno));
        }
    }
}


string CTar::x_ToFilesystemPath(const string& name) const
{
    string path(CDirEntry::IsAbsolutePath(name)  ||  m_BaseDir.empty()
                ? name : CDirEntry::ConcatPath(m_BaseDir, name));
    return CDirEntry::NormalizePath(path);
}


string CTar::x_ToArchiveName(const string& path) const
{
    // NB: Path assumed to have been normalized
    string retval = CDirEntry::AddTrailingPathSeparator(path);

#ifdef NCBI_OS_MSWIN
    // Convert to Unix format with forward slashes
    NStr::ReplaceInPlace(retval, "\\", "/");
    const NStr::ECase how = NStr::eNocase;
#else
    const NStr::ECase how = NStr::eCase;
#endif // NCBI_OS_MSWIN

    bool absolute;
    // Remove leading base dir from the path
    if (!m_BaseDir.empty()  &&  NStr::StartsWith(retval, m_BaseDir, how)) {
        if (retval.length() > m_BaseDir.length()) {
            retval.erase(0, m_BaseDir.length()/*separator too*/);
        } else {
            retval.assign(".", 1);
        }
        absolute = false;
    } else {
        absolute = CDirEntry::IsAbsolutePath(retval);
    }

    SIZE_TYPE pos = 0;

#ifdef NCBI_OS_MSWIN
    _ASSERT(!retval.empty());
    // Remove disk name if present
    if (isalpha((unsigned char) retval[0])  &&  retval.find(":") == 1) {
        pos = 2;
    }
#endif // NCBI_OS_MSWIN

    // Remove any leading and trailing slashes
    while (pos < retval.length()  &&  retval[pos] == '/') {
        pos++;
    }
    if (pos) {
        retval.erase(0, pos);
    }
    pos = retval.length();
    while (pos > 0  &&  retval[pos - 1] == '/') {
        --pos;
    }
    if (pos < retval.length()) {
        retval.erase(pos);
    }

    // Check on '..'
    if (retval == ".."  ||  NStr::StartsWith(retval, "../")  ||
        NStr::EndsWith(retval, "/..")  ||  retval.find("/../") != NPOS) {
        TAR_THROW(this, eBadName,
                  "Name '" + retval + "' embeds parent directory ('..')");
    }

    if (absolute) {
        retval.insert((SIZE_TYPE) 0, 1, '/');
    }
    return retval;
}


auto_ptr<CTar::TEntries> CTar::x_Append(const string&   name,
                                        const TEntries* toc)
{
    auto_ptr<TEntries> entries(new TEntries);

    const EFollowLinks follow_links = (m_Flags & fFollowLinks ?
                                       eFollowLinks : eIgnoreLinks);

    // Compose entry name for relative names
    string path = x_ToFilesystemPath(name);

    // Get direntry information
    CDirEntry entry(path);
    CDirEntry::SStat st;
    if (!entry.Stat(&st, follow_links)) {
        int x_errno = errno;
        TAR_THROW(this, eOpen,
                  "Cannot get status of '" + path + '\''+ s_OSReason(x_errno));
    }
    CDirEntry::EType type = CDirEntry::GetType(st.orig);

    // Create the entry info
    m_Current = CTarEntryInfo(m_StreamPos);
    m_Current.m_Name = x_ToArchiveName(path);
    if (type == CDirEntry::eDir  &&  m_Current.m_Name != "/") {
        m_Current.m_Name += '/';
    }
    if (m_Current.GetName().empty()) {
        TAR_THROW(this, eBadName,
                  "Empty entry name not allowed");
    }
    m_Current.m_Type = CTarEntryInfo::EType(type);
    if (m_Current.GetType() == CTarEntryInfo::eSymLink) {
        _ASSERT(!follow_links);
        m_Current.m_LinkName = entry.LookupLink();
        if (m_Current.GetLinkName().empty()) {
            TAR_THROW(this, eBadName,
                      "Empty link name not allowed");
        }
    }
    entry.GetOwner(&m_Current.m_UserName, &m_Current.m_GroupName);
    m_Current.m_Stat = st.orig;
    // Fixup for mode bits
    m_Current.m_Stat.st_mode = (mode_t) s_ModeToTar(st.orig.st_mode);

    // Check if we need to update this entry in the archive
    bool update = true;
    if (toc) {
        bool found = false;

        if (type != CDirEntry::eUnknown) {
            // Start searching from the end of the list, to find
            // the most recent entry (if any) first
            string temp;
            REVERSE_ITERATE(TEntries, i, *toc) {
                if (!temp.empty()) {
                    if (i->GetType() == CTarEntryInfo::eHardLink
                        ||  temp != x_ToFilesystemPath(i->GetName())) {
                        continue;
                    }
                } else if (path == x_ToFilesystemPath(i->GetName())) {
                    found = true;
                    if (i->GetType() == CTarEntryInfo::eHardLink) {
                        temp = x_ToFilesystemPath(i->GetLinkName());
                        continue;
                    }
                } else {
                    continue;
                }
                if (m_Current.GetType() != i->GetType()) {
                    if (m_Flags & fEqualTypes) {
                        goto out;
                    }
                } else if (m_Current.GetType() == CTarEntryInfo::eSymLink
                           &&  m_Current.GetLinkName() == i->GetLinkName()) {
                    goto out;
                }
                if (m_Current.GetModificationTime() <=
                    i->GetModificationTime()) {
                    update = false; // same(or older), no update
                }
                break;
            }
        }

        if (!update  ||  (!found  &&  (m_Flags & (fUpdate & ~fOverwrite)))) {
            if (type != CDirEntry::eDir  &&  type != CDirEntry::eUnknown) {
                goto out;
            }
            // Directories always get recursive treatment later
            update = false;
        }
    }

    // Append entry
    switch (type) {
    case CDirEntry::eFile:
        _ASSERT(update);
        x_AppendFile(path);
        entries->push_back(m_Current);
        break;

    case CDirEntry::eBlockSpecial:
    case CDirEntry::eCharSpecial:
    case CDirEntry::ePipe:
    case CDirEntry::eLink:
        _ASSERT(update);
        /*FALLTHRU*/
    case CDirEntry::eDir:
        if (update) {
            m_Current.m_Stat.st_size = 0;
            x_WriteEntryInfo(path);
            entries->push_back(m_Current);
        }
        if (type == CDirEntry::eDir) {
            // Append/Update all files from that directory
            CDir::TEntries dir = CDir(path).GetEntries("*",
                                                       CDir::eIgnoreRecursive);
            ITERATE(CDir::TEntries, i, dir) {
                auto_ptr<TEntries> e = x_Append((*i)->GetPath(), toc);
                entries->splice(entries->end(), *e.get());
            }
        }
        break;

    case CDirEntry::eDoor:
    case CDirEntry::eSocket:
        // Tar does not have any provisions to store this kind of entry
        TAR_POST(3, Warning,
                 "Skipping non-archiveable "
                 + string(type == CDirEntry::eSocket ? "socket" : "door")
                 + " entry '" + path + '\'');
        break;

    case CDirEntry::eUnknown:
        TAR_THROW(this, eBadName,
                  "Unable to handle '" + path + '\'');
        /*FALLTHRU*/

    default:
        TAR_POST(14, Warning,
                 "Skipping unsupported source '" + path
                 + "' w/type #" + NStr::IntToString(int(type)));
        break;
    }

 out:
    return entries;
}


// Regular files only!
void CTar::x_AppendFile(const string& file)
{
    // FIXME:  Switch to CFileIO eventually to avoid ifstream
    // obscurity w.r.t. errors, an extra layer of buffering etc.
    CNcbiIfstream ifs;

    _ASSERT(m_Current.GetType() == CTarEntryInfo::eFile);

    // Open file
    ifs.open(file.c_str(), IOS_BASE::binary | IOS_BASE::in);
    if (!ifs) {
        int x_errno = errno;
        TAR_THROW(this, eOpen,
                  "Cannot open file '" + file + '\'' + s_OSReason(x_errno));
    }

    // Write file header
    x_WriteEntryInfo(file);

    Uint8 size = m_Current.GetSize();
    while (size) {
        // Write file contents
        size_t avail = m_BufferSize - m_BufferPos;
        if (avail > size) {
            avail = (size_t) size;
        }
        // Read file
        if (!ifs.read(m_Buffer + m_BufferPos, (streamsize) avail)) {
            int x_errno = errno;
            TAR_THROW(this, eRead,
                      "Error reading file '" +file+ '\''+ s_OSReason(x_errno));
        }
        avail = (size_t) ifs.gcount();
        // Write buffer to the archive
        x_WriteArchive(avail);
        size -= avail;
    }

    // Write zeros to get the written size a multiple of BLOCK_SIZE
    size_t zero = ALIGN_SIZE(m_BufferPos) - m_BufferPos;
    memset(m_Buffer + m_BufferPos, 0, zero);
    x_WriteArchive(zero);
    _ASSERT(!OFFSET_OF(m_BufferPos));
}


class CTarReader : public IReader
{
public:
    CTarReader(CTar* tar, EOwnership own = eNoOwnership)
        : m_Read(0), m_Eof(false), m_Bad(false), m_Own(own), m_Tar(tar)
    { }
    ~CTarReader();

    virtual ERW_Result Read(void* buf, size_t count, size_t* bytes_read = 0);
    virtual ERW_Result PendingCount(size_t* count);

private:
    Uint8      m_Read;
    bool       m_Eof;
    bool       m_Bad;
    EOwnership m_Own;
    CTar*      m_Tar;
};


CTarReader::~CTarReader()
{
    if (m_Own) {
        delete m_Tar;
    }
}


ERW_Result CTarReader::Read(void* buf, size_t count, size_t* bytes_read)
{
    if (m_Bad  ||  !count) {
        if (bytes_read) {
            *bytes_read = 0;
        }
        return m_Bad ? eRW_Error
            : (m_Read < m_Tar->m_Current.GetSize()  ||  !m_Eof) ? eRW_Success
            : eRW_Eof;
    }

    size_t read;
    if (m_Read >= m_Tar->m_Current.GetSize()) {
        m_Eof = true;
        read = 0;
    } else {
        size_t left = (size_t)(m_Tar->m_Current.GetSize() - m_Read);
        if (count > left) {
            count = left;
        }

        size_t off = (size_t) OFFSET_OF(m_Read);
        if (off) {
            read = BLOCK_SIZE - off;
            if (m_Tar->m_BufferPos) {
                off += m_Tar->m_BufferPos  - BLOCK_SIZE;
            } else {
                off += m_Tar->m_BufferSize - BLOCK_SIZE;
            }
            if (read > count) {
                read = count;
            }
            memcpy(buf, m_Tar->m_Buffer + off, read);
            m_Read += read;
            count  -= read;
            if (!count) {
                goto out;
            }
            buf = (char*) buf + read;
        } else {
            read = 0;
        }

        off = m_Tar->m_BufferPos;
        if (m_Tar->x_ReadArchive(count)) {
            m_Tar->m_StreamPos += ALIGN_SIZE(count);
            memcpy(buf, m_Tar->m_Buffer + off, count);
            m_Read += count;
            read   += count;
        } else {
            m_Bad = true;
            TAR_THROW(m_Tar, eRead,
                      "Read error while streaming");
        }
    }

 out:
    if (bytes_read) {
        *bytes_read = read;
    }
    return m_Bad ? eRW_Error : m_Eof ? eRW_Eof : eRW_Success;
}


ERW_Result CTarReader::PendingCount(size_t* count)
{
    if (m_Bad) {
        return eRW_Error;
    }
    Uint8 left = m_Tar->m_Current.GetSize() - m_Read;
    if (!left  &&  m_Eof) {
        return eRW_Eof;
    }
    *count = (size_t) left;
    return eRW_Success;
}


IReader* CTar::Extract(istream& is, const string& name, CTar::TFlags flags)
{
    auto_ptr<CTar> tar(new CTar(is, 1/*blocking factor*/));
    tar->SetFlags(flags);

    auto_ptr<CMaskFileName> mask(new CMaskFileName);
    mask->Add(name);
    tar->SetMask(mask.get(), eTakeOwnership);
    mask.release();
   
    auto_ptr<TEntries> temp = tar->x_Open(eInternal);
    _ASSERT(temp.get()  &&  temp->size() < 2);
    if (temp->size() < 1) {
        return 0;
    }

    _ASSERT(tar->m_Current == *temp->begin());
    CTarEntryInfo::EType type = tar->m_Current.GetType();
    if (type != CTarEntryInfo::eFile  &&
        (type != CTarEntryInfo::eUnknown  ||  (flags & fSkipUnsupported))) {
        return 0;
    }

    IReader* ir = new CTarReader(tar.get(), eTakeOwnership);
    tar.release();
    return ir;
}


IReader* CTar::GetNextEntryData(void)
{
    return
        m_Current.GetType() != CTarEntryInfo::eFile ? 0 : new CTarReader(this);
}


END_NCBI_SCOPE
