#ifndef CORELIB___RWSTREAMBUF__HPP
#define CORELIB___RWSTREAMBUF__HPP

/*  $Id: rwstreambuf.hpp 195071 2010-06-18 16:37:44Z kazimird $
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
 *   Reader-writer based stream buffer
 *
 */

/// @file rwstreambuf.hpp
/// Reader-writer based stream buffer
/// @sa IReader, IWriter, IReaderWriter

#include <corelib/ncbistre.hpp>
#include <corelib/reader_writer.hpp>

#ifdef NCBI_COMPILER_MIPSPRO
#  define CRWStreambufBase CMIPSPRO_ReadsomeTolerantStreambuf
#else
#  define CRWStreambufBase CNcbiStreambuf
#  ifdef NCBI_COMPILER_MSVC
#    pragma warning(push)
#    pragma warning(disable:4996)
#  endif //NCBI_COMPILER_MSVC
#endif //NCBI_COMPILER_MIPSPRO


BEGIN_NCBI_SCOPE


/// Reader-writer based stream buffer

class NCBI_XNCBI_EXPORT CRWStreambuf : public CRWStreambufBase
{
public:
    /// Which of the objects (passed in the constructor) should be
    /// deleted on this object's destruction.
    /// NOTE:  if the reader and writer are in fact the same object,
    ///        it will _not_ be deleted twice.
    enum EFlags {
        fOwnReader      = 1 << 1,    // own the underlying reader
        fOwnWriter      = 1 << 2,    // own the underlying writer
        fOwnAll         = fOwnReader + fOwnWriter,
        fLogExceptions  = 1 << 8,
        fLeakExceptions = 1 << 9
    };
    typedef int TFlags;             // bitwise OR of EFlags


    CRWStreambuf(IReaderWriter* rw       = 0,
                 streamsize     buf_size = 0,
                 CT_CHAR_TYPE*  buf      = 0,
                 TFlags         flags    = 0);

    CRWStreambuf(IReader*       r,
                 IWriter*       w,
                 streamsize     buf_size = 0,
                 CT_CHAR_TYPE*  buf      = 0,
                 TFlags         flags    = 0);

    virtual ~CRWStreambuf();

protected:
    virtual CT_INT_TYPE overflow(CT_INT_TYPE c);
    virtual streamsize  xsputn(const CT_CHAR_TYPE* buf, streamsize n);

    virtual CT_INT_TYPE underflow(void);
    virtual streamsize  xsgetn(CT_CHAR_TYPE* s, streamsize n);
    virtual streamsize  showmanyc(void);

    virtual int         sync(void);

    /// Note: setbuf(0, 0) has no effect
    virtual CNcbiStreambuf* setbuf(CT_CHAR_TYPE* buf, streamsize buf_size);

    // only seekoff(0, IOS_BASE::cur, *) is permitted
    virtual CT_POS_TYPE seekoff(CT_OFF_TYPE off, IOS_BASE::seekdir whence,
                                IOS_BASE::openmode which =
                                IOS_BASE::in | IOS_BASE::out);
protected:
    CT_POS_TYPE    x_GetGPos(void)
    { return x_GPos - (CT_OFF_TYPE)(gptr() ? egptr() - gptr() : 0); }
    CT_POS_TYPE    x_GetPPos(void)
    { return x_PPos + (CT_OFF_TYPE)(pptr() ? pbase() - pptr() : 0); }

protected:
    TFlags         m_Flags;

    IReader*       m_Reader;
    IWriter*       m_Writer;

    streamsize     m_BufSize;
    CT_CHAR_TYPE*  m_ReadBuf;
    CT_CHAR_TYPE*  m_WriteBuf;

    CT_CHAR_TYPE*  m_pBuf;
    CT_CHAR_TYPE   x_Buf;

    CT_POS_TYPE    x_GPos;      //< get position [for istream.tellg()]
    CT_POS_TYPE    x_PPos;      //< put position [for ostream.tellp()]

    bool           x_Err;       //< whether there was a write error
    CT_POS_TYPE    x_ErrPos;    //< position of write error (if x_Error)
};


END_NCBI_SCOPE


#ifdef NCBI_COMPILER_MSVC
#  pragma warning(pop)
#endif //NCBI_COMPILER_MSVC

#endif /* CORELIB___RWSTREAMBUF__HPP */
