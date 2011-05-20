/*  $Id: reader_base.hpp 167543 2009-08-04 11:54:16Z ludwigf $
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
 * Author: Frank Ludwig
 *
 * File Description:
 *   Basic reader interface
 *
 */

#ifndef OBJTOOLS_READERS___READERBASE__HPP
#define OBJTOOLS_READERS___READERBASE__HPP

#include <corelib/ncbistd.hpp>
#include <objects/seq/Seq_annot.hpp>
#include <util/format_guess.hpp>

BEGIN_NCBI_SCOPE

BEGIN_objects_SCOPE // namespace ncbi::objects::

class IErrorContainer;
class CObjReaderLineException;

//  ----------------------------------------------------------------------------
class NCBI_XOBJREAD_EXPORT CReaderBase
//  ----------------------------------------------------------------------------
{
public:
    enum ObjectType {
        OT_UNKNOWN,
        OT_SEQANNOT,
        OT_SEQENTRY
    };
    
public:

    virtual ~CReaderBase() {}

    //
    //  Class interface:
    //
    static CReaderBase* GetReader(
        CFormatGuess::EFormat,
        int =0 );

    //
    //  Object interface:
    //
    virtual CRef< CSerialObject >
    ReadObject(
        CNcbiIstream&,
        IErrorContainer* =0 );
                
    virtual CRef< CSerialObject >
    ReadObject(
        ILineReader&,
        IErrorContainer* =0 ) =0;
                
    virtual CRef< CSeq_annot >
    ReadSeqAnnot(
        CNcbiIstream&,
        IErrorContainer* =0 );
                
    virtual CRef< CSeq_annot >
    ReadSeqAnnot(
        ILineReader&,
        IErrorContainer* =0 );
                
    virtual CRef< CSeq_entry >
    ReadSeqEntry(
        CNcbiIstream&,
        IErrorContainer* =0 );
                
    virtual CRef< CSeq_entry >
    ReadSeqEntry(
        ILineReader&,
        IErrorContainer* =0 );
                
    //
    //  Class helper functions:
    //
protected:
    virtual bool x_ParseBrowserLine(
        const string&,
        CRef<CSeq_annot>& );
        
    virtual bool x_ParseTrackLine(
        const string&,
        CRef<CSeq_annot>& );
        
    virtual void x_GetTrackValues(
        const string&,
        map<string,string>& );
        
    virtual void x_SetBrowserRegion(
        const string&,
        CAnnot_descr& );

    virtual void x_SetTrackData(
        CRef<CSeq_annot>&,
        CRef<CUser_object>&,
        const string&,
        const string& );
                
    virtual void x_AddConversionInfo(
        CRef< CSeq_annot >&,
        IErrorContainer* );
                    
    virtual void x_AddConversionInfo(
        CRef< CSeq_entry >&,
        IErrorContainer* );
                    
    static bool SplitLines( 
        const char* pcBuffer, 
        size_t uSize,
        list<string>& lines );

    void
    ProcessError(
        CObjReaderLineException&,
        IErrorContainer* );
        
    //
    //  Data:
    //
    unsigned int m_uLineNumber;
    int m_iFlags;
};

END_objects_SCOPE
END_NCBI_SCOPE

#endif // OBJTOOLS_READERS___READERBASE__HPP
