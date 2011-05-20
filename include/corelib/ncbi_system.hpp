#ifndef NCBI_SYSTEM__HPP
#define NCBI_SYSTEM__HPP

/*  $Id: ncbi_system.hpp 220394 2011-01-19 23:05:29Z rafanovi $
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
 * Author:  Vladimir Ivanov
 *
 * File Description: System functions
 *
 */


#include <corelib/ncbitime.hpp>

#if defined(NCBI_OS_MSWIN)
#  include <corelib/ncbi_os_mswin.hpp>
#endif


BEGIN_NCBI_SCOPE


/////////////////////////////////////////////////////////////////////////////
///
/// Process limits
///
/// If, during the program execution, the process exceed any from limits
/// (see ELimitsExitCodeMemory) then:
///   1) Dump info about current program's state to log-stream.
///      - if defined print handler "handler", then it will be used.
///      - if defined "parameter", it will be passed into print handler;
///   2) Terminate the program.
/// One joint print handler for all limit types is used.


/// Code for program's exit handler.
enum ELimitsExitCode {
    eLEC_None,    ///< Normal exit.
    eLEC_Memory,  ///< Memory limit.
    eLEC_Cpu      ///< CPU usage limit.
};

/// Type of parameter for print handler.
typedef void* TLimitsPrintParameter;

/// Type of handler for printing a dump information after generating
/// any limitation event.
typedef void (*TLimitsPrintHandler)(ELimitsExitCode, size_t, CTime&, TLimitsPrintParameter);


/// [UNIX only]  Set memory limit.
///
/// Set the limit for the size of dynamic memory (heap) allocated
/// by the process. 
///
/// @note 
///   The implementation of malloc() can be different. Some systems use 
///   sbrk()-based implementation, other use mmap() system call to allocate
///   memory, ignoring data segment. In the second case SetHeapLimit() 
///   don't work at all. Usually don't know about how exactly malloc()
///   is implemented. We added another function - SetMemoryLimit(), that
///   supports mmap()-based memory allocation, please use it instead.
/// 
/// @param max_heap_size
///   The maximal amount of dynamic memory can be allocated by the process.
///   (including heap)
///   The 0 value lift off the heap restrictions.
/// @param handler
///   Pointer to a print handler used for dump output.
///   Use default handler if passed as NULL.
/// @param parameter
///   Parameter carried into the print handler. Can be passed as NULL.
/// @return 
///   Completion status.
/// @sa SetMemoryLimit
/// @deprecated
NCBI_DEPRECATED
NCBI_XNCBI_EXPORT
extern bool SetHeapLimit(size_t max_size, 
                         TLimitsPrintHandler handler = 0, 
                         TLimitsPrintParameter parameter = 0);


/// [UNIX only]  Set memory limit.
///
/// Set the limit for the size of used memory allocated by the process.
/// 
/// @param max_size
///   The maximal amount of memory in bytes that can be allocated by
///   the process. Use the same limits for process's data segment
///   (including heap) and virtual memory (address space).
///   On 32-bit systems limit is at most 2 GiB, or this resource is unlimited.
///   The 0 value lift off the heap restrictions.
/// @param handler
///   Pointer to a print handler used for dump output in the case of reaching
///   memory limit. Use default handler if passed as NULL.
/// @param parameter
///   Parameter carried into the print handler. Can be passed as NULL.
///   Usefull if singular handler is used for setting some limits.
///   See also SetCpuTimeLimit().
/// @return 
///   Completion status.
/// @note
///   Setting limits may not work on some systems, depends on OS, compilation
///   options and etc. Some systems enforce memory limits, other didn't.
/// @sa SetCpuTimeLimit, TLimitsPrintHandler
NCBI_XNCBI_EXPORT
extern bool SetMemoryLimit(size_t max_size, 
                           TLimitsPrintHandler handler = 0, 
                           TLimitsPrintParameter parameter = 0);


/// [UNIX only]  Set CPU usage limit.
///
/// Set the limit for the CPU time that can be consumed by current process.
/// 
/// @param max_cpu_time
///   The maximal amount of seconds of CPU time can be consumed by the process.
///   The 0 value lifts off the CPU time restrictions if allowed to do so.
/// @param handler
///   Pointer to a print handler used for dump output in the case of reaching
///   CPU usage limit. Use default handler if passed as NULL.
/// @param parameter
///   Parameter carried into the print handler. Can be passed as NULL.
/// @terminate_time
///   The time in seconds that the process will have to terminate itself after
///   receiving a signal about exceeding CPU usage limit. After that it can
///   be killed by OS.
/// @return 
///   Completion status.
/// @note
///   Setting a low CPU time limit cannot be generally undone to a value
///   higher than "max_cpu_time + terminate_time" at a later time.
/// @sa SetMemoryLimit, TLimitsPrintHandler
NCBI_XNCBI_EXPORT
extern bool SetCpuTimeLimit(size_t                max_cpu_time,
                            TLimitsPrintHandler   handler = 0, 
                            TLimitsPrintParameter parameter = 0,
                            size_t                terminate_time = 5);


/////////////////////////////////////////////////////////////////////////////
///
/// System information
///

/// [UNIX & Windows]  Return number of active CPUs (never less than 1).
NCBI_XNCBI_EXPORT
extern unsigned int GetCpuCount(void);


/// [UNIX & Windows]
/// Return granularity with which virtual memory is allocated.
/// Return 0 if cannot determine it on current platform or if an error occurs.
NCBI_XNCBI_EXPORT
extern unsigned long GetVirtualMemoryPageSize(void);

/// [UNIX & Windows]
/// Return the amount of physical memory available on the machine.
/// Return 0 if cannot determine it on current platform or if an error occurs.
NCBI_XNCBI_EXPORT
extern Uint8 GetPhysicalMemorySize(void);

/// [UNIX & Windows] Return current memory usage, in bytes.
/// NULL arguments will not be filled in.
/// Returns true if able to determine memory usage, and false otherwise.
NCBI_XNCBI_EXPORT
extern bool GetMemoryUsage(size_t* total, size_t* resident, size_t* shared);

/////////////////////////////////////////////////////////////////////////////
///
/// Sleep
///
/// Suspend execution for a time.
///
/// Sleep for at least the specified number of microsec/millisec/seconds.
/// Time slice restrictions are imposed by platform/OS.
/// On UNIX the sleep can be interrupted by a signal.
/// Sleep*Sec(0) have no effect (but may cause context switches).
///
/// [UNIX & Windows]

NCBI_XNCBI_EXPORT
extern void SleepSec(unsigned long sec,
                     EInterruptOnSignal onsignal = eRestartOnSignal);

NCBI_XNCBI_EXPORT
extern void SleepMilliSec(unsigned long ml_sec,
                          EInterruptOnSignal onsignal = eRestartOnSignal);

NCBI_XNCBI_EXPORT
extern void SleepMicroSec(unsigned long mc_sec,
                          EInterruptOnSignal onsignal = eRestartOnSignal);



/////////////////////////////////////////////////////////////////////////////
///
/// Suppress Diagnostic Popup Messages
///

/// Suppress modes
enum ESuppressSystemMessageBox {
    fSuppress_System    = (1<<0),     ///< System errors
    fSuppress_Runtime   = (1<<1),     ///< Runtime library
    fSuppress_Debug     = (1<<2),     ///< Debug library
    fSuppress_Exception = (1<<3),     ///< Unhandled exceptions
    fSuppress_All       = fSuppress_System | fSuppress_Runtime | 
                          fSuppress_Debug  | fSuppress_Exception,
    fSuppress_Default   = fSuppress_All
};
/// Binary OR of "ESuppressSystemMessageBox"
typedef int TSuppressSystemMessageBox;  

/// Suppress popup messages on execution errors.
///
/// NOTE: MS Windows-specific.
/// Suppresses all error message boxes in both runtime and in debug libraries,
/// as well as all General Protection Fault messages.
NCBI_XNCBI_EXPORT
extern void SuppressSystemMessageBox(TSuppressSystemMessageBox mode = 
                                     fSuppress_Default);

/// Prevent run of SuppressSystemMessageBox().
///
/// NOTE: MS Windows-specific.
/// If this function is called, all following calls of
/// SuppressSystemMessageBox() will be ignored. If SuppressSystemMessageBox()
/// was executed before, that this function print out a critical error message.
/// For example can be used in CGI applications where SuppressSystemMessageBox
/// always calls in the CCgiApplication constructor.
/// 
NCBI_XNCBI_EXPORT
extern void DisableSuppressSystemMessageBox();


/// Check if system message box has been suppressed for debug library.
///
/// NOTE: MS Windows-specific.
NCBI_XNCBI_EXPORT
extern bool IsSuppressedDebugSystemMessageBox();

END_NCBI_SCOPE

#endif  /* NCBI_SYSTEM__HPP */
