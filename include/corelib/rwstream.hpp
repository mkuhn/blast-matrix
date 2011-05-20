#ifndef CORELIB___RWSTREAM__HPP
#define CORELIB___RWSTREAM__HPP

/*  $Id: rwstream.hpp 257881 2011-03-16 15:14:53Z rafanovi $
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
 * Authors:  Anton Lavrentiev
 *
 * File Description:
 *   Reader-writer based streams
 *
 */

/// @file rwstream.hpp
/// Reader-writer based streams
/// @sa IReader, IWriter, IReaderWriter, CRWStreambuf


#include <corelib/ncbimisc.hpp>
#include <corelib/impl/rwstreambuf.hpp>


BEGIN_NCBI_SCOPE


/** @addtogroup Stream
 *
 * @{
 */


/// Note about the "buf_size" parameter for the streams in this API.
///
/// CRWStream implementation is targeted at minimizing in-memory data
/// copy operations assosiated with I/O for intermediate buffering.
/// For that, the following policies apply:
///
/// 1.  No read operation from the output device shall be requested less than
///     the specified stream buffer size.  In cases when the user code requests
///     more than the internal stream buffer size, the request may be passed
///     through to the input device to read directly into the user provided
///     buffer.  In that case, the read request can only be larger than the
///     the size of the internal stream buffer.
///
/// 2.  Write operations from an input device are done in full buffers, unless:
///  a. An incoming user write request is larger than the internal buffer,
///     then the contents of the internal buffer gets flushed first
///     (which may comprise of fewer than the buffer size bytes) followed
///     by the direct write request of the user's block (larger than the
///     internal stream buffer);
///  b. Flushing of an internal buffer (including 2a above) resulted in a
///     short write on the device (fewer bytes actually written), then the
///     successive write attempt may contain fewer bytes than the size of
///     the internal stream buffer (namely, the remainder of what has been
///     left behind by the preceding write attempt).
///
/// However, any portable implementation should *not* rely on how data chunks
/// are being flushed or requested by the stream implementations.  If further
/// factoring into blocks (e.g. specifically-sized) is necessary for an I/O
/// device to operate properly, that should be implemented at the level of
/// respective IReader/IWriter API explicitly.


/// Reader-based stream; @sa IReader
///
/// @param buf_size
///     specifies the number bytes for internal I/O buffer, entirely used
///     for reading by the underlying stream buffer CRWStreambuf;
///     0 causes to create the buffer of some default size.
///
/// @param buf
///     may specify the buffer location (if 0, an internal storage gets
///     allocated and later freed upon stream destruction).
///
/// @param flags
///     controls whether IReader is destroyed upon stream destruction,
///     and whether exceptions cause logging (or caught silently).
///
/// Special case of "buf_size" == 1 and "buf" == 0 creates unbuffered stream.
///
/// @sa IWStream, IRWStream

class NCBI_XNCBI_EXPORT CRStream : public CNcbiIstream
{
public:
    CRStream(IReader*             r,
             streamsize           buf_size = 0,
             CT_CHAR_TYPE*        buf      = 0,
             CRWStreambuf::TFlags flags    = 0) :
        CNcbiIstream(0), m_Sb(r, 0, buf_size, buf, flags)
    {
        init(&m_Sb);
    }

private:
    CRWStreambuf m_Sb;
};


/// Writer-based stream; @sa IWriter
///
/// @param buf_size
///     specifies the number bytes for internal I/O buffer, entirely used
///     for writing by the underlying stream buffer CRWStreambuf;
///     0 causes to create the buffer of some default size.
///
/// @param buf
///     may specify the buffer location (if 0, an internal storage gets
///     allocated and later freed upon stream destruction).
///
/// @param flags
///     controls whether IWriter is destroyed upon stream destruction,
///     and whether exceptions cause logging (or caught silently).
///
/// Special case of "buf_size" == 1 and "buf" == 0 creates unbuffered stream.
///
/// @sa IRStream, IRWStream

class NCBI_XNCBI_EXPORT CWStream : public CNcbiOstream
{
public:
    CWStream(IWriter*             w,
             streamsize           buf_size = 0,
             CT_CHAR_TYPE*        buf      = 0,
             CRWStreambuf::TFlags flags    = 0) :
        CNcbiOstream(0), m_Sb(0, w, buf_size, buf, flags)
    {
        init(&m_Sb);
    }

private:
    CRWStreambuf m_Sb;
};


/// Reader-writer based stream; @sa IReaderWriter
///
/// @param buf_size
///     specifies the number bytes for internal I/O buffer,
///     with half used for reading and the other half for writing
///     by underlying stream buffer CRWStreambuf;
///     0 causes to create the buffer of some default size.
///
/// @param buf
///     may specify the buffer location (if 0, an internal storage gets
///     allocated and later freed upon stream destruction).
///
/// @param flags
///     controls whether IReader is destroyed upon stream destruction,
///     and whether exceptions cause logging (or caught silently).
///
/// Special case of "buf_size" == 1 and "buf" == 0 creates unbuffered stream.
///
/// @sa IRStream, IWStream

class NCBI_XNCBI_EXPORT CRWStream : public CNcbiIostream
{
public:
    CRWStream(IReaderWriter*       rw,
              streamsize           buf_size = 0,
              CT_CHAR_TYPE*        buf      = 0,
              CRWStreambuf::TFlags flags    = 0)
        : CNcbiIostream(0), m_Sb(rw, buf_size, buf, flags)
    {
        init(&m_Sb);
    }

private:
    CRWStreambuf m_Sb;
};


/// istream-based IReader
class NCBI_XNCBI_EXPORT CStreamReader : public IReader
{
public:
    CStreamReader(CNcbiIstream& is, EOwnership own = eNoOwnership)
        : m_Stream(&is, own)
    {
    }
    ~CStreamReader();

    virtual ERW_Result Read(void* buf, size_t count, size_t* bytes_read = 0);
    virtual ERW_Result PendingCount(size_t* count);

private:
    AutoPtr<CNcbiIstream> m_Stream;

private: // prevent copy
    CStreamReader(const CStreamReader&);
    void operator=(const CStreamReader&);
};


/* @} */


END_NCBI_SCOPE

#endif /* CORELIB___RWSTREAM__HPP */
