/* $Id: ncbi_conn_streambuf.cpp 257812 2011-03-16 14:56:04Z rafanovi $
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
 * Author:  Anton Lavrentiev, Denis Vakatov
 *
 * File Description:
 *   CONN-based C++ stream buffer
 *
 */

#include <ncbi_pch.hpp>
#include "ncbi_conn_streambuf.hpp"
#include <corelib/ncbidbg.hpp>
#include <connect/ncbi_conn_exception.hpp>
#include <connect/error_codes.hpp>
#include <memory>
#include <stdlib.h>  // free()

#define NCBI_USE_ERRCODE_X   Connect_Stream


BEGIN_NCBI_SCOPE


CConn_Streambuf::CConn_Streambuf(CONNECTOR       connector,
                                 const STimeout* timeout,
                                 streamsize      buf_size,
                                 bool            tie,
                                 CT_CHAR_TYPE*   ptr,
                                 size_t          size)
    : m_Conn(0), m_WriteBuf(0), m_BufSize(buf_size ? buf_size : 1),
      m_Status(eIO_Unknown), m_Tie(tie), m_Close(true), m_CbValid(false),
      x_GPos((CT_OFF_TYPE)(ptr ? size : 0)), x_PPos((CT_OFF_TYPE) size)
{
    if ( !connector ) {
        ERR_POST_X(2, x_Message("CConn_Streambuf(): NULL connector"));
        return;
    }
    if ((m_Status = CONN_Create(connector, &m_Conn)) != eIO_Success) {
        ERR_POST_X(3, x_Message("CConn_Streambuf(): CONN_Create() failed"));
        if (connector->destroy)
            connector->destroy(connector);
        return;
    }
    _ASSERT(m_Conn != 0);
    x_Init(timeout, buf_size, ptr, size);
}


CConn_Streambuf::CConn_Streambuf(CONN            conn,
                                 bool            close,
                                 const STimeout* timeout,
                                 streamsize      buf_size,
                                 bool            tie,
                                 CT_CHAR_TYPE*   ptr,
                                 size_t          size)
    : m_Conn(conn), m_WriteBuf(0), m_BufSize(buf_size ? buf_size : 1),
      m_Status(eIO_Unknown), m_Tie(tie), m_Close(close), m_CbValid(false),
      x_GPos((CT_OFF_TYPE)(ptr ? size : 0)), x_PPos((CT_OFF_TYPE) size)
{
    if ( !m_Conn ) {
        ERR_POST_X(1, x_Message("CConn_Streambuf(): NULL connection"));
        return;
    }
    x_Init(timeout, buf_size, ptr, size);
}


void CConn_Streambuf::x_Init(const STimeout* timeout, streamsize buf_size,
                             CT_CHAR_TYPE* ptr, size_t size)
{
    if (timeout != kDefaultTimeout) {
        CONN_SetTimeout(m_Conn, eIO_Open,      timeout);
        CONN_SetTimeout(m_Conn, eIO_ReadWrite, timeout);
        CONN_SetTimeout(m_Conn, eIO_Close,     timeout);
    }

    m_WriteBuf = buf_size ? new CT_CHAR_TYPE[m_BufSize << 1] : 0;
    m_ReadBuf  = buf_size ? m_WriteBuf + m_BufSize           : &x_Buf;

    setp(m_WriteBuf,    m_WriteBuf + buf_size);  // Put area (if any)
    if ( ptr )
        setg(ptr,       ptr,       ptr + size);  // Initial get area
    else
        setg(m_ReadBuf, m_ReadBuf, m_ReadBuf);   // Empty get area

    SCONN_Callback cb;
    cb.func = x_OnClose; /* NCBI_FAKE_WARNING: WorkShop */
    cb.data = this;
    CONN_SetCallback(m_Conn, eCONN_OnClose, &cb, &m_Cb);
    m_CbValid = true;
    m_Status = eIO_Success;
}


void CConn_Streambuf::x_Close(bool close)
{
    if ( !m_Conn )
        return;

    // Flush only if data pending
    if (pbase()  &&  pptr() > pbase())
        sync();

    setg(0, 0, 0);
    setp(0, 0);

    CONN c = m_Conn;
    m_Conn = 0;  // NB: no re-entry

    if ( close ) {
        // here: not called from close callback
        if ( m_CbValid ) {
            SCONN_Callback cb;
            CONN_SetCallback(c, eCONN_OnClose, &m_Cb, &cb);
            if ((void*) cb.func != (void*) x_OnClose  ||  cb.data != this)
                CONN_SetCallback(c, eCONN_OnClose, &cb, 0);
        }

        if (m_Close  &&  (m_Status = CONN_Close(c)) != eIO_Success)
            _TRACE(x_Message("x_Close(): CONN_Close() failed"));
    } else if (m_CbValid  &&  m_Cb.func)
        m_Cb.func(c, eCONN_OnClose, m_Cb.data);
}


void CConn_Streambuf::x_OnClose(CONN           _DEBUG_ARG(conn),
                                ECONN_Callback _DEBUG_ARG(type),
                                void*          data)
{
    CConn_Streambuf* sb = static_cast<CConn_Streambuf*>(data);

    _ASSERT(type == eCONN_OnClose  &&  sb  &&  conn);
    _ASSERT(!sb->m_Conn  ||  sb->m_Conn == conn);
    sb->x_Close(false);
}


string CConn_Streambuf::x_Message(const char* msg)
{
    const char* type = m_Conn ? CONN_GetType(m_Conn)     : 0;
    char* descr      = m_Conn ? CONN_Description(m_Conn) : 0;
    string result("CConn_Streambuf::");
    result += msg;
    result += " (";
    result += type ? type : "UNKNOWN";
    if (descr) {
        result += "; ";
        result += descr;
        free(descr);
    }
    result += "): ";
    result += IO_StatusStr(m_Status);
    return result;
}


CT_INT_TYPE CConn_Streambuf::overflow(CT_INT_TYPE c)
{
    if ( !m_Conn )
        return CT_EOF;

    if ( pbase() ) {
        // send buffer
        size_t n_towrite = pptr() - pbase();
        if ( n_towrite ) {
            size_t n_written;
            m_Status = CONN_Write(m_Conn, pbase(),
                                  n_towrite, &n_written, eIO_WritePlain);
            if ( !n_written ) {
                _ASSERT(m_Status != eIO_Success);
                ERR_POST_X(4, x_Message("overflow(): CONN_Write() failed"));
                return CT_EOF;
            }
            // update buffer content (get rid of the data just sent)
            memmove(pbase(), pbase() + n_written, n_towrite - n_written);
            x_PPos += (CT_OFF_TYPE) n_written;
            pbump(-int(n_written));
        }

        // store char
        if ( !CT_EQ_INT_TYPE(c, CT_EOF) )
            return sputc(CT_TO_CHAR_TYPE(c));
    } else if ( !CT_EQ_INT_TYPE(c, CT_EOF) ) {
        // send char
        size_t n_written;
        CT_CHAR_TYPE b = CT_TO_CHAR_TYPE(c);
        m_Status = CONN_Write(m_Conn, &b, 1, &n_written, eIO_WritePlain);
        if ( !n_written ) {
            _ASSERT(m_Status != eIO_Success);
            ERR_POST_X(5, x_Message("overflow(): CONN_Write(1) failed"));
        }
        x_PPos += (CT_OFF_TYPE) n_written;
        return n_written == 1 ? c : CT_EOF;
    }

    _ASSERT(CT_EQ_INT_TYPE(c, CT_EOF));
    if ((m_Status = CONN_Flush(m_Conn)) != eIO_Success) {
        ERR_POST_X(9, x_Message("overflow(): CONN_Flush() failed"));
        return CT_EOF;
    }
    return CT_NOT_EOF(CT_EOF);
}


streamsize CConn_Streambuf::xsputn(const CT_CHAR_TYPE* buf, streamsize m)
{
    if ( !m_Conn )
        return 0;

    if (m <= 0)
        return 0;

    size_t n = (size_t) m;
    size_t n_written = 0;

    m_Status = eIO_Success;
    for (;;) {
        size_t x_written;

        if ( pbase() ) {
            if (pbase() + n < epptr()) {
                // Would entirely fit into the buffer not causing an overflow
                x_written = (size_t)(epptr() - pptr());
                if (x_written > n)
                    x_written = n;
                if ( x_written ) {
                    memcpy(pptr(), buf, x_written);
                    pbump(int(x_written));
                    n_written += x_written;
                    n         -= x_written;
                    if (!n)
                        return n_written;
                    buf       += x_written;
                }
            }
            if (m_Status != eIO_Success)
                break;

            size_t x_towrite = (size_t)(pptr() - pbase());
            if ( x_towrite ) {
                m_Status = CONN_Write(m_Conn, pbase(),
                                      x_towrite, &x_written, eIO_WritePlain);
                if ( !x_written ) {
                    _ASSERT(m_Status != eIO_Success);
                    ERR_POST_X(6, x_Message("xsputn():"
                                            " CONN_Write(Plain) failed"));
                    break;
                }
                memmove(pbase(), pbase() + x_written, x_towrite - x_written);
                x_PPos += (CT_OFF_TYPE) x_written;
                pbump(-int(x_written));
                continue;
            }
        }

        _ASSERT(n  &&  m_Status == eIO_Success);
        m_Status = CONN_Write(m_Conn, buf, n, &x_written, eIO_WritePersist);
        if ( !x_written ) {
            _ASSERT(m_Status != eIO_Success);
            ERR_POST_X(7, x_Message("xsputn():"
                                    " CONN_Write(Persist) failed"));
            if ( !pbase() )
                return (streamsize) n_written;
            break;
        }
        x_PPos    += (CT_OFF_TYPE) x_written;
        n_written += x_written;
        n         -= x_written;
        if (!n  ||  !pbase())
            return (streamsize) n_written;
        buf       += x_written;
        if (m_Status != eIO_Success)
            break;
    }

    _ASSERT(n  &&  pbase());
    if (pptr() < epptr()) {
        size_t x_written = (size_t)(epptr() - pptr());
        if (x_written > n)
            x_written = n;
        memcpy(pptr(), buf, x_written);
        n_written += x_written;
        pbump(int(x_written));
    }

    return (streamsize) n_written;
}


CT_INT_TYPE CConn_Streambuf::underflow(void)
{
    if ( !m_Conn )
        return CT_EOF;

    // flush output buffer, if tied up to it
    if (m_Tie  &&  sync() != 0)
        return CT_EOF;

    _ASSERT(!gptr()  ||  gptr() >= egptr());

#ifdef NCBI_COMPILER_MIPSPRO
    if (m_MIPSPRO_ReadsomeGptrSetLevel  &&  m_MIPSPRO_ReadsomeGptr != gptr())
        return CT_EOF;
    m_MIPSPRO_ReadsomeGptr = (CT_CHAR_TYPE*)(-1L);
#endif /*NCBI_COMPILER_MIPSPRO*/

    // read from connection
    size_t n_read;
    m_Status = CONN_Read(m_Conn, m_ReadBuf,
                         m_BufSize, &n_read, eIO_ReadPlain);
    if ( !n_read ) {
        _ASSERT(m_Status != eIO_Success);
        if (m_Status != eIO_Closed)
            ERR_POST_X(8, x_Message("underflow(): CONN_Read() failed"));
        return CT_EOF;
    }

    // update input buffer with the data just read
    x_GPos += (CT_OFF_TYPE) n_read;
    setg(m_ReadBuf, m_ReadBuf, m_ReadBuf + n_read);

    return CT_TO_INT_TYPE(*m_ReadBuf);
}


streamsize CConn_Streambuf::xsgetn(CT_CHAR_TYPE* buf, streamsize m)
{
    if ( !m_Conn )
        return 0;

    // flush output buffer, if tied up to it
    if (m_Tie  &&  sync() != 0)
        return 0;

    if (m <= 0)
        return 0;
    size_t n = (size_t) m;

    // first, read from the memory buffer
    size_t n_read = gptr() ? (size_t)(egptr() - gptr()) : 0;
    if (n_read > n)
        n_read = n;
    memcpy(buf, gptr(), n_read);
    gbump(int(n_read));
    buf += n_read;
    n   -= n_read;

    while ( n ) {
        // next, read from the connection
        size_t       x_read = n < (size_t) m_BufSize ? m_BufSize : n;
        CT_CHAR_TYPE* x_buf = n < (size_t) m_BufSize ? m_ReadBuf : buf;
        m_Status = CONN_Read(m_Conn, x_buf, x_read, &x_read, eIO_ReadPlain);
        if ( !x_read ) {
            _ASSERT(m_Status != eIO_Success);
            if (m_Status != eIO_Closed)
                ERR_POST_X(10, x_Message("xsgetn(): CONN_Read() failed"));
            break;
        }
        x_GPos += (CT_OFF_TYPE) x_read;
        // satisfy "usual backup condition", see standard: 27.5.2.4.3.13
        if (x_buf == m_ReadBuf) {
            size_t xx_read = x_read;
            if (x_read > n)
                x_read = n;
            memcpy(buf, m_ReadBuf, x_read);
            setg(m_ReadBuf, m_ReadBuf + x_read, m_ReadBuf + xx_read);
        } else {
            _ASSERT(x_read <= n);
            size_t xx_read = x_read > (size_t) m_BufSize ? m_BufSize : x_read;
            memcpy(m_ReadBuf, buf + x_read - xx_read, xx_read);
            setg(m_ReadBuf, m_ReadBuf + xx_read, m_ReadBuf + xx_read);
        }
        n_read += x_read;
        if (m_Status != eIO_Success)
            break;
        buf    += x_read;
        n      -= x_read;
    }

    return (streamsize) n_read;
}


streamsize CConn_Streambuf::showmanyc(void)
{
    static const STimeout kZeroTmo = {0, 0};

    if ( !m_Conn )
        return -1;

    // flush output buffer, if tied up to it
    if (m_Tie  &&  sync() != 0)
        return -1;

    _ASSERT(!gptr()  ||  gptr() >= egptr());

    const STimeout* tmo;
    const STimeout* timeout = CONN_GetTimeout(m_Conn, eIO_Read);
    if (timeout == kDefaultTimeout) {
        // HACK * HACK * HACK
        tmo = ((SMetaConnector*) m_Conn)->default_timeout;
    } else
        tmo = timeout;

    size_t x_read;
    if (m_BufSize > 1) {
        if ( !tmo )
            _VERIFY(CONN_SetTimeout(m_Conn, eIO_Read, &kZeroTmo)==eIO_Success);
        m_Status = CONN_Read(m_Conn, m_ReadBuf + 1, m_BufSize - 1,
                             &x_read, eIO_ReadPlain);
        if ( !tmo )
            _VERIFY(CONN_SetTimeout(m_Conn, eIO_Read, timeout)  ==eIO_Success);
        _ASSERT(x_read > 0  ||  m_Status != eIO_Success);
    } else {
        x_read = 0;
        m_Status = CONN_Wait(m_Conn, eIO_Read, tmo ? tmo : &kZeroTmo);
    }

    if ( !x_read ) {
        switch (m_Status) {
        case eIO_Success:
            _ASSERT(m_BufSize <= 1);
            return  1;  // can read at least 1 byte
        case eIO_Timeout:
            if (!tmo  ||  !(tmo->sec | tmo->usec))
                break;
            /*FALLTHRU*/
        case eIO_Closed:
            return -1;  // EOF
        default:
            break;
        }
        return      0;  // no data available immediately
    }

    _ASSERT(m_BufSize > 1);
    if ( gptr() )
        m_ReadBuf[0] = gptr()[-1];
    setg(m_ReadBuf, m_ReadBuf + 1, m_ReadBuf + 1 + x_read);
    x_GPos += x_read;
    return x_read;
}


int CConn_Streambuf::sync(void)
{
    if ( !m_Conn )
        return -1;

    do {
        if (CT_EQ_INT_TYPE(overflow(CT_EOF), CT_EOF))
            return -1;
    } while (pbase()  &&  pptr() > pbase());
    return 0;
}


CNcbiStreambuf* CConn_Streambuf::setbuf(CT_CHAR_TYPE* /*buf*/,
                                        streamsize    /*buf_size*/)
{
    NCBI_THROW(CConnException, eConn, "CConn_Streambuf::setbuf() not allowed");
    /*NOTREACHED*/
    return this; /*dummy for compiler*/ /* NCBI_FAKE_WARNING */
}


CT_POS_TYPE CConn_Streambuf::seekoff(CT_OFF_TYPE        off,
                                     IOS_BASE::seekdir  whence,
                                     IOS_BASE::openmode which)
{
    if (m_Conn  &&  off == 0  &&  whence == IOS_BASE::cur) {
        // tellp()/tellg() support only
        switch (which) {
        case IOS_BASE::out:
            return x_PPos + (CT_OFF_TYPE)(pptr() ? pptr() - pbase() : 0);
        case IOS_BASE::in:
            return x_GPos - (CT_OFF_TYPE)(gptr() ? egptr() - gptr() : 0);
        default:
            break;
        }
    }
    return (CT_POS_TYPE)((CT_OFF_TYPE)(-1));
}


const char* CConnException::GetErrCodeString(void) const
{
    switch (GetErrCode()) {
    case eConn:    return "eConn";
    default:       return CException::GetErrCodeString();
    }
}


END_NCBI_SCOPE
