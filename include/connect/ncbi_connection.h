#ifndef CONNECT___NCBI_CONNECTION__H
#define CONNECT___NCBI_CONNECTION__H

/* $Id: ncbi_connection.h 257880 2011-03-16 15:14:36Z rafanovi $
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
 * Author:  Denis Vakatov
 *
 * File Description:
 *   Generic API to open and handle connection to an abstract I/O service.
 *   Several methods can be used to establish the connection, and each of them
 *   yields in a simple handle(of type "CONN") that contains a handle(of type
 *   "CONNECTOR") to a data and methods implementing the generic connection I/O
 *   operations. E.g. this API can be used to:
 *     1) connect using HTTPD-based dispatcher (e.g. to NCBI services);
 *     2) hit a CGI script;
 *     3) connect to a bare socket at some "host:port";
 *     4) whatever else can fit this paradigm -- see the SConnectorTag-related
 *        structures;  e.g. it could be a plain file I/O or even a memory area.
 *
 *  See in "ncbi_connector.h" for the detailed specification of the underlying
 *  connector("CONNECTOR", "SConnectorTag") methods and data structures.
 *
 */

#include <connect/ncbi_connector.h>


/** @addtogroup Connectors
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


struct SConnectionTag;
typedef struct SConnectionTag* CONN;      /* connection handle */


/* Compose all data necessary to establish a new connection
 * (merely bind it to the specified connector). Unsuccessful completion
 * sets conn to 0, and leaves connector intact (can be used again).
 * NOTE1:  The real connection will not be established right away. Instead,
 *         it will be established at the moment of the first call to one of
 *         "Flush", "Wait", "Write", or "Read" methods.
 * NOTE2:  "Connection establishment" at this level of abstraction may differ
 *         from actual link establishment at the underlying connector's level.
 * NOTE3:  Initial timeout values are set to kDefaultTimeout, meaning
 *         that connector-specific timeouts are in force for the connection.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Create
(CONNECTOR connector,  /* [in]  connector                        */
 CONN*     conn        /* [out] handle of the created connection */
 );


/* Reinit using new "connector".
 * If "conn" is already opened then close the current connection at first,
 * even if "connector" is just the same as the current connector.
 * If "connector" is NULL then close and destroy the incumbent,
 * and leave connection empty (effective way to destroy connector(s)).
 * NOTE:  Although it closes the previous connection immediately, however it
 *        does not open the new connection right away:  see notes on "Create".
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_ReInit
(CONN      conn,      /* [in] connection handle */
 CONNECTOR connector  /* [in] new connector     */
 );


/* Get verbal representation of connection type as a character string.
 * Note that the returned value is only valid until the next
 * I/O operation in the connection.  Return value NULL denotes
 * unknown connection type.
 */
extern NCBI_XCONNECT_EXPORT const char* CONN_GetType
(CONN conn  /* [in]  connection handle */ 
 );


/* Get read (event == eIO_Read) or write (event == eIO_Write)
 * position within the connection.
 * Positions are advanced from 0 on, and only concerning I/O that has
 * caused calling to the actual connector's "read" (i.e. pushbacks
 * never considered, and peeks -- not always) and "write" methods.
 * Special case:  eIO_Open as "event" causes to clear both positions
 * with 0, and to return 0.
 */
extern NCBI_XCONNECT_EXPORT TNCBI_BigCount CONN_GetPosition
(CONN      conn,  /* [in]  connection handle */ 
 EIO_Event event  /* [in]  see description   */
 );


/* Return human-readable description of the connection as a character
 * '\0'-terminated string.  The string is not guaranteed to have any
 * particular format and is intended solely for something like
 * logging and debugging.  Return NULL if the connection cannot
 * provide any description information (or if it is in a bad state).
 * Application program must call free() to deallocate space occupied
 * by the returned string when the description is no longer needed.
 */
extern NCBI_XCONNECT_EXPORT char* CONN_Description
(CONN conn  /* [in]  connection handle */
 );


/* Specify timeout for the connection I/O, including "Connect" (aka "Open")
 * and "Close".  May be called at any time during the connection lifetime.
 * NOTE1:  if "new_timeout" is NULL then set the timeout to be infinite.
 * NOTE2:  if "new_timeout" is kDefaultTimeout then an underlying,
 *         connector-specific value is used (this is the default).
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_SetTimeout
(CONN            conn,        /* [in] connection handle */
 EIO_Event       event,       /* [in] I/O direction     */
 const STimeout* new_timeout  /* [in] new timeout       */
 );


/* Retrieve current timeout (return NULL if it is infinite).
 * The returned pointer is guaranteed to point to a valid timeout structure,
 * or to be either NULL or kDefaultTimeout until next "SetTimeout"
 * or "Close" method's call.
 */
extern NCBI_XCONNECT_EXPORT const STimeout* CONN_GetTimeout
(CONN      conn,  /* [in] connection handle                  */
 EIO_Event event  /* [in] I/O direction, not "eIO_ReadWrite" */
 );


/* Block on the connection until it becomes available for either read or
 * write (dep. on "event"), until timeout expires, or until any error.
 * NOTE:  "timeout" can also be one of two special values:
 *         NULL (means infinite), kDefaultTimeout (connector-defined).
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Wait
(CONN            conn,    /* [in] connection handle                  */
 EIO_Event       event,   /* [in] can be eIO_Read or eIO_Write only! */
 const STimeout* timeout  /* [in] the maximal wait time              */
 );


/* Write up to "size" bytes from the buffer "buf" to the connection.
 * Return the number of actually written bytes in "*n_written".
 * It may not return "eIO_Success" if no data at all can be written before
 * write timeout expired or an error occurred.
 * Parameter "how" modifies the write behavior:
 * eIO_WritePlain   -- return immediately after having written as many
 *                     as 1 byte of data, or if an error has occurred;
 * eIO_WritePersist -- return only after having written all of the data
 *                     from "buf", or if an error has occurred.
 * NOTE:  See CONN_SetTimeout() how to set write timeout.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Write
(CONN            conn,      /* [in]  connection handle                     */ 
 const void*     buf,       /* [in]  pointer to the data buffer to write   */ 
 size_t          size,      /* [in]  # of bytes to write                   */ 
 size_t*         n_written, /* [out, non-NULL] # of actually written bytes */
 EIO_WriteMethod how        /* [in]  eIO_WritePlain or eIO_WritePersist    */
 );


/* Push back "size" bytes from the buffer "buf" into connection.
 * Return eIO_Success on success, other code on error.
 * NOTE1:  Data pushed back may not necessarily be the same as obtained
 *         from the connection before.
 * NOTE2:  Upon following read operation, the pushed back data are
 *         taken out first.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_PushBack
(CONN        conn,  /* [in]  connection handle                     */
 const void* buf,   /* [in]  pointer to the data being pushed back */
 size_t      size   /* [in]  # of bytes to push back               */
 );


/* Explicitly flush connection from any pending data written by "CONN_Write()".
 * NOTE1:  CONN_Flush() effectively opens connection (if it wasn't open yet).
 * NOTE2:  Connection considered open if underlying connector's "Open" method
 *         has successfully executed; actual data link may not yet exist.
 * NOTE3:  CONN_Read() always calls CONN_Flush() before proceeding;
 *         so does CONN_Close() but only if connection is was open before.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Flush
(CONN        conn   /* [in] connection handle                      */
 );


/* Read up to "size" bytes from connection to the buffer pointed to by "buf".
 * Return the number of actually read bytes in "*n_read".
 * May not return eIO_Success if no data at all can be read before
 * read timeout expired or an error occurred.
 * Parameter "how" modifies the read behavior:
 *   eIO_ReadPlain   -- return immediately after having read as many as
 *                      1 byte from connection, or if an error has occurred;
 *   eIO_ReadPeek    -- eIO_ReadPlain but don't discard read data from CONN;
 *   eIO_ReadPersist -- return only after having filled full "buf" with data
 *                      (exactly "size" bytes), or if an error has occurred.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Read
(CONN           conn,   /* [in]  connection handle                  */
 void*          buf,    /* [out] memory buffer to read to           */
 size_t         size,   /* [in]  max. # of bytes to read            */
 size_t*        n_read, /* [out, non-NULL] # of actually read bytes */
 EIO_ReadMethod how     /* [in]  read/peek | persist                */
 );


/* Read up to "size" bytes from connection into the string buffer pointed
 * to by "line".  Stop reading if either '\n' or an error is encountered.
 * Replace '\n' with '\0'.  Upon return "*n_read" contains the number
 * of characters written to "line", not including the terminating '\0'.
 * If not enough space provided in "line" to accomodate the '\0'-terminated
 * line, then all "size" bytes are used and "*n_read" equals "size" on return.
 * This is the only case when "line" will not be '\0'-terminated.
 * Return code advises the caller whether another read can be attempted:
 *   eIO_Success -- read completed successfully, keep reading;
 *   other code  -- an error occurred, and further attempt may fail.
 *
 * This call utilizes eIO_Read timeout as set by CONN_SetTimeout().
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_ReadLine
(CONN    conn,
 char*   line,
 size_t  size,
 size_t* n_read
 );


/* Obtain status of the last IO operation. This is NOT a completion
 * code of the last CONN-call, but rather a status from the lower level
 * connector's layer.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Status
(CONN      conn,   /* [in]  connection handle       */
 EIO_Event dir     /* [in] = {eIO_Read | eIO_Write} */
 );


/* Cancel the connection's I/O ability.
 * This is *not* connection closure, but any data extraction or
 * insertion (Read/Write) will be effectively rejected after this call
 * (and eIO_Interrupt will result, same for CONN_Status()).
 * CONN_Close() is still required to release internal connection
 * structures.
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Cancel
(CONN conn  /* [in] connection handle */
 );


/* Close the connection, destroy relevant internal data.
 * NOTE:  whatever error code is returned, the connection handle "conn"
 *        will become invalid (so, you should not use it anymore).
 */
extern NCBI_XCONNECT_EXPORT EIO_Status CONN_Close
(CONN conn  /* [in] connection handle */
 );


/* Set user callback function to be called upon an event specified by the
 * callback type.  Note that the callback function is always called prior
 * to the event to happen, e.g. the eCONN_OnClose callback is called when
 * the connection is about to close, but have not yet been closed.
 * The callback function is supplied with 3 arguments: connection handle,
 * type of event, and the user data (specified when the callback was set).
 * CONN_SetCallback() stores previous callback in "old_cb" (if it is not NULL).
 * The callbacks are called only once (they get reset each time prior to
 * been actually called), so the code that wants to get callbacks repeatedly
 * must reinstate them as necessary with CONN_SetCallback() calls
 * (e.g. from inside the callbacks themselves).
 */
typedef enum {
    eCONN_OnClose  = 0,  /* NB: connection has been flushed prior to the call*/
    eCONN_OnRead   = 1,  /* Read from connector is about to occur            */
    eCONN_OnWrite  = 2,  /* Write to connector is about to occur             */
    eCONN_OnCancel = 3   /* CONN_Cancel() is about to take effect            */
} ECONN_Callback;
#define CONN_N_CALLBACKS 4

typedef void (*FConnCallback)(CONN conn, ECONN_Callback type, void* data);

typedef struct {
    FConnCallback func;  /* Function to call on the event                */
    void*         data;  /* Data to pass to the callback as its last arg */
} SCONN_Callback;

extern NCBI_XCONNECT_EXPORT EIO_Status CONN_SetCallback
(CONN                  conn,    /* [in]  connection to set callback for      */
 ECONN_Callback        type,    /* [in]  callback type                       */
 const SCONN_Callback* new_cb,  /* [in]  callback to set (may be 0 to reset) */
 SCONN_Callback*       old_cb   /* [out] to save old callback at (may be 0)  */
);


#ifdef __cplusplus
}  /* extern "C" */
#endif


/* @} */

#endif /* CONNECT___NCBI_CONNECTION__H */
