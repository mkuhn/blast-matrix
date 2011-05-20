/*  $Id: seq_entry_ci.cpp 103491 2007-05-04 17:18:18Z kazimird $
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
* Author: Aleksey Grichenko, Eugene Vasilchenko
*
* File Description:
*    Handle to Seq-entry object
*
*/

#include <ncbi_pch.hpp>
#include <objmgr/seq_entry_ci.hpp>

#include <objmgr/impl/seq_entry_info.hpp>
#include <objmgr/impl/scope_impl.hpp>
#include <objmgr/impl/bioseq_set_info.hpp>

#include <objects/seqset/Bioseq_set.hpp>


BEGIN_NCBI_SCOPE
BEGIN_SCOPE(objects)


/////////////////////////////////////////////////////////////////////////////
// CSeq_entry_CI
/////////////////////////////////////////////////////////////////////////////


CSeq_entry_CI::CSeq_entry_CI(const CSeq_entry_Handle& entry)
{
    x_Initialize(entry.GetSet());
}


CSeq_entry_CI::CSeq_entry_CI(const CBioseq_set_Handle& seqset)
{
    x_Initialize(seqset);
}


void CSeq_entry_CI::x_Initialize(const CBioseq_set_Handle& seqset)
{
    if ( seqset ) {
        m_Parent = seqset;
        m_Iterator = seqset.x_GetInfo().GetSeq_set().begin();
        x_SetCurrentEntry();
    }
}


void CSeq_entry_CI::x_SetCurrentEntry(void)
{
    if ( m_Parent && m_Iterator != m_Parent.x_GetInfo().GetSeq_set().end() ) {
        m_Current = CSeq_entry_Handle(**m_Iterator,
                                      m_Parent.GetTSE_Handle());
    }
    else {
        m_Current.Reset();
    }
}


CSeq_entry_CI& CSeq_entry_CI::operator ++(void)
{
    if ( *this ) {
        ++m_Iterator;
        x_SetCurrentEntry();
    }
    return *this;
}


/////////////////////////////////////////////////////////////////////////////
// CSeq_entry_I
/////////////////////////////////////////////////////////////////////////////


CSeq_entry_I::CSeq_entry_I(const CSeq_entry_EditHandle& entry)
{
    x_Initialize(entry.SetSet());
}


CSeq_entry_I::CSeq_entry_I(const CBioseq_set_EditHandle& seqset)
{
    x_Initialize(seqset);
}


void CSeq_entry_I::x_Initialize(const CBioseq_set_EditHandle& seqset)
{
    if ( seqset ) {
        m_Parent = seqset;
        m_Iterator = seqset.x_GetInfo().SetSeq_set().begin();
        x_SetCurrentEntry();
    }
}


void CSeq_entry_I::x_SetCurrentEntry(void)
{
    if ( m_Parent && m_Iterator != m_Parent.x_GetInfo().SetSeq_set().end() ) {
        m_Current = CSeq_entry_EditHandle(**m_Iterator,
                                          m_Parent.GetTSE_Handle());
    }
    else {
        m_Current.Reset();
    }
}


CSeq_entry_I& CSeq_entry_I::operator ++(void)
{
    if ( *this ) {
        ++m_Iterator;
        x_SetCurrentEntry();
    }
    return *this;
}


END_SCOPE(objects)
END_NCBI_SCOPE
