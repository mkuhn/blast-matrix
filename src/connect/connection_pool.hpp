#ifndef CONNECT___CONNECTION_POOL__HPP
#define CONNECT___CONNECTION_POOL__HPP

/* $Id: connection_pool.hpp 194770 2010-06-17 10:35:19Z ivanovp $
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
 * Authors:  Aaron Ucko, Victor Joukov
 *
 */

/// @file connection_pool.hpp
/// Internal header for threaded server connection pools.


#include <connect/impl/server_connection.hpp>


/** @addtogroup ThreadedServer
 *
 * @{
 */


BEGIN_NCBI_SCOPE


class CServer_ControlConnection : public CSocket,
                                  public IServer_ConnectionBase
{
public:
    virtual CStdRequest* CreateRequest(EServIO_Event event,
                                       CServer_ConnectionPool& connPool,
                                       const STimeout* timeout,
                                       int request_id);
};


class CServer_ConnectionPool
{
public:
    CServer_ConnectionPool(unsigned max_connections);
    ~CServer_ConnectionPool();

    typedef IServer_ConnectionBase TConnBase;
    typedef CServer_Connection     TConnection;
    typedef CServer_Listener       TListener;

    enum EConnType {
        eInactiveSocket,
        eActiveSocket,
        eListener,
        ePreDeferredSocket,
        eDeferredSocket,
        ePreClosedSocket,
        eClosedSocket
    };

    void SetMaxConnections(unsigned max_connections) {
        m_MaxConnections = max_connections;
    }

    bool Add(TConnBase* conn, EConnType type);
    void Remove(TConnBase* conn);
    void PingControlConnection(void);

    /// Guard connection from out-of-order packet processing by
    /// pulling eActiveSocket's from poll vector
    /// Resets the expiration time as a bonus.
    void SetConnType(TConnBase* conn, EConnType type);

    /// Close connection as if it was initiated by server (not by client).
    void CloseConnection(TConnBase* conn);

    /// Clean up inactive connections which are no longer open or
    /// which have been idle for too long.
    void Clean(vector<IServer_ConnectionBase*>& revived_conns);

    /// Erase all connections
    void Erase(void);

    bool GetPollAndTimerVec(vector<CSocketAPI::SPoll>& polls,
        vector<IServer_ConnectionBase*>& timer_requests,
        STimeout* timer_timeout) const;

    void StartListening(void);
    void StopListening(void);

    struct SPerConnInfo {
        void UpdateExpiration(const TConnBase* conn);

        CTime     expiration;
        EConnType type;
    };

private:
    typedef map<TConnBase*, SPerConnInfo> TData;
    TData                                 m_Data;
    mutable CMutex                        m_Mutex;
    unsigned int                          m_MaxConnections;
    CSocket                               m_ControlSocket;
    mutable CServer_ControlConnection     m_ControlSocketForPoll;
};


END_NCBI_SCOPE


/* @} */

#endif  /* CONNECT___CONNECTION_POOL__HPP */
