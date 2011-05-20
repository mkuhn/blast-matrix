/*  $Id: exception.cpp 103491 2007-05-04 17:18:18Z kazimird $
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
* Author: Eugene Vasilchenko
*
* File Description:
*   Standard serialization exceptions
*/

#include <ncbi_pch.hpp>
#include <serial/exception.hpp>

BEGIN_NCBI_SCOPE

void CSerialException::AddFrameInfo(string frame_info)
{
    m_FrameStack = frame_info + m_FrameStack;
}

void CSerialException::ReportExtra(ostream& out) const
{
    if ( !m_FrameStack.empty() ) {
        out << " at " << m_FrameStack;
    }
}

const char* CInvalidChoiceSelection::GetErrCodeString(void) const
{
    switch ( GetErrCode() ) {
    case eFail:  return "eFail";
    default:     return CException::GetErrCodeString();
    }
}

const char* CInvalidChoiceSelection::GetName(
    size_t index, const char* const names[], size_t namesCount)
{
    if ( index > namesCount )
        return "?unknown?";
    return names[index];
    
}

CInvalidChoiceSelection::CInvalidChoiceSelection(
    const CDiagCompileInfo& diag_info,
    size_t currentIndex, size_t mustBeIndex,
    const char* const names[], size_t namesCount, 
    EDiagSev severity)
        : CSerialException(diag_info, 0,
          (CSerialException::EErrCode) CException::eInvalid,"")
{
    x_Init(diag_info,
           string("Invalid choice selection: ")+
           GetName(currentIndex, names, namesCount)+". "
           "Expected: "+
           GetName(mustBeIndex, names, namesCount),0, severity);
    x_InitErrCode((CException::EErrCode)(CInvalidChoiceSelection::eFail));
}

CInvalidChoiceSelection::CInvalidChoiceSelection(
    const char* file, int line,
    size_t currentIndex, size_t mustBeIndex,
    const char* const names[], size_t namesCount,
    EDiagSev severity)
        : CSerialException(CDiagCompileInfo(file, line), 0,
          (CSerialException::EErrCode) CException::eInvalid,"")
{
    x_Init(CDiagCompileInfo(file, line),
           string("Invalid choice selection: ")+
           GetName(currentIndex, names, namesCount)+". "
           "Expected: "+
           GetName(mustBeIndex, names, namesCount),0, severity);
    x_InitErrCode((CException::EErrCode)(CInvalidChoiceSelection::eFail));
}

CInvalidChoiceSelection::CInvalidChoiceSelection(
    size_t currentIndex, size_t mustBeIndex,
    const char* const names[], size_t namesCount,
    EDiagSev severity)
        : CSerialException(CDiagCompileInfo("unknown", 0), 0,
          (CSerialException::EErrCode) CException::eInvalid,"")
{
    x_Init(CDiagCompileInfo("unknown", 0),
           string("Invalid choice selection: ")+
           GetName(currentIndex, names, namesCount)+". "
           "Expected: "+
           GetName(mustBeIndex, names, namesCount),0, severity);
    x_InitErrCode((CException::EErrCode)(CInvalidChoiceSelection::eFail));
}

CInvalidChoiceSelection::CInvalidChoiceSelection(
    const CInvalidChoiceSelection& other)
    : CSerialException(other)
{
    x_Assign(other);
}

CInvalidChoiceSelection::~CInvalidChoiceSelection(void) throw()
{
}

const char* CInvalidChoiceSelection::GetType(void) const
{
    return "CInvalidChoiceSelection";
}

CInvalidChoiceSelection::TErrCode CInvalidChoiceSelection::GetErrCode(void) const
{
    return typeid(*this) == typeid(CInvalidChoiceSelection) ?
        (CInvalidChoiceSelection::TErrCode) x_GetErrCode() :
        (CInvalidChoiceSelection::TErrCode) CException::eInvalid;
}

CInvalidChoiceSelection::CInvalidChoiceSelection(void)
{
}

const CException* CInvalidChoiceSelection::x_Clone(void) const
{
    return new CInvalidChoiceSelection(*this);
}

const char* CSerialException::GetErrCodeString(void) const
{
    switch ( GetErrCode() ) {
    case eNotImplemented: return "eNotImplemented";
    case eEOF:            return "eEOF";
    case eIoError:        return "eIoError";
    case eFormatError:    return "eFormatError";
    case eOverflow:       return "eOverflow";
    case eInvalidData:    return "eInvalidData";
    case eIllegalCall:    return "eIllegalCall";
    case eFail:           return "eFail";
    case eNotOpen:        return "eNotOpen";
    case eMissingValue:   return "eMissingValue";
    default:              return CException::GetErrCodeString();
    }
}

const char* CUnassignedMember::GetErrCodeString(void) const
{
#if 0
    switch ( GetErrCode() ) {
        case eGet:            return "eGet";
        case eWrite:          return "eWrite";
        case eUnknownMember:  return "eUnknownMember";
        default:              return CException::GetErrCodeString();
        }
#else
        // At least with ICC 9.0 on 64-bit Linux in Debug and MT mode
    // there is an apparent bug that causes the above "switch" based
    // variant of this function to misbehave and crash with SEGV...
    TErrCode e = GetErrCode();
    if (e == eGet)           {return "eGet";}
    if (e == eWrite)         {return "eWrite";}
    if (e == eUnknownMember) {return "eUnknownMember";}
    return CException::GetErrCodeString();
#endif
}


END_NCBI_SCOPE
