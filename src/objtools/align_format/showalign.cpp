/*  $Id: showalign.cpp 227410 2011-01-27 21:51:57Z rafanovi $
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
 *  and reliability of the software and data, the NLM and thesubset U.S.
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
 * Author:  Jian Ye
 *
 * File Description:
 *   Sequence alignment display
 *
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = "$Id: showalign.cpp 227410 2011-01-27 21:51:57Z rafanovi $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <ncbi_pch.hpp>

#include <objtools/align_format/showalign.hpp>

#include <corelib/ncbiexpt.hpp>
#include <corelib/ncbiutil.hpp>
#include <corelib/ncbistre.hpp>
#include <corelib/ncbireg.hpp>

#include <util/range.hpp>
#include <util/md5.hpp>

#include <objects/general/Object_id.hpp>
#include <objects/general/User_object.hpp>
#include <objects/general/User_field.hpp>
#include <objects/general/Dbtag.hpp>

#include <serial/iterator.hpp>
#include <serial/objistr.hpp>
#include <serial/objostr.hpp>
#include <serial/serial.hpp>
#include <serial/objostrasnb.hpp> 
#include <serial/objistrasnb.hpp> 
#include <connect/ncbi_conn_stream.hpp>

#include <objmgr/scope.hpp>
#include <objmgr/feat_ci.hpp>
#include <objtools/data_loaders/genbank/gbloader.hpp>

#include <objmgr/util/sequence.hpp>
#include <objmgr/util/feature.hpp>

#include <objects/seqfeat/SeqFeatData.hpp>
#include <objects/seqfeat/Cdregion.hpp>
#include <objects/seqfeat/Genetic_code.hpp>
#include <objects/seq/Seq_descr.hpp>
#include <objects/seq/Seqdesc.hpp>
#include <objects/seq/Bioseq.hpp>

#include <objects/seqset/Seq_entry.hpp>

#include <objects/seqloc/Seq_id.hpp>
#include <objects/seqloc/Seq_interval.hpp>

#include <objects/seqalign/Seq_align_set.hpp>
#include <objects/seqalign/Score.hpp>
#include <objects/seqalign/Std_seg.hpp>
#include <objects/seqalign/Dense_diag.hpp>

#include <objtools/alnmgr/alnmix.hpp>
#include <objtools/alnmgr/alnvec.hpp>

#include <objects/blastdb/Blast_def_line.hpp>
#include <objects/blastdb/Blast_def_line_set.hpp>
#include <objtools/blast/seqdb_reader/linkoutdb.hpp>

#include <stdio.h>
#include <util/tables/raw_scoremat.h>
#include <objtools/readers/getfeature.hpp>
#include <html/htmlhelper.hpp>


BEGIN_NCBI_SCOPE
USING_SCOPE(objects);
USING_SCOPE(sequence);
BEGIN_SCOPE(align_format)

static const char k_IdentityChar = '.';
static const int k_NumFrame = 6;
static const string k_FrameConversion[k_NumFrame] = {"+1", "+2", "+3", "-1",
                                                     "-2", "-3"};
static const int k_GetSubseqThreshhold = 10000;

///threshhold to color mismatch. 98 means 98% 
static const int k_ColorMismatchIdentity = 0; 
static const int k_GetDynamicFeatureSeqLength = 200000;
static const string k_DumpGnlUrl = "/blast/dumpgnl.cgi";
static const int k_FeatureIdLen = 16;
const string color[]={"#000000", "#808080", "#FF0000"};
const string k_ColorRed = "#FF0000";
const string k_ColorPink = "#F805F5";

static const char k_IntronChar = '~';
static const int k_IdStartMargin = 2;
static const int k_SeqStopMargin = 2;
static const int k_StartSequenceMargin = 2;

static const string k_UncheckabeCheckbox = "<input type=\"checkbox\" \
name=\"getSeqMaster\" value=\"\" onClick=\"uncheckable('getSeqAlignment%d',\
 'getSeqMaster')\">";

static const string k_Checkbox = "<input type=\"checkbox\" \
name=\"getSeqGi\" value=\"%s\" onClick=\"synchronizeCheck(this.value, \
'getSeqAlignment%d', 'getSeqGi', this.checked)\">";

static const string k_CheckboxEx = "<input type=\"checkbox\" name=\"getSeqGi\" value=\"%s\" \
checked=\"checked\" onClick=\"synchAl(this);\">";
#ifdef USE_ORG_IMPL
static string k_GetSeqSubmitForm[] = {"<FORM  method=\"post\" \
action=\"http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?SUBMIT=y\" \
name=\"%s%d\"><input type=button value=\"Get selected sequences\" \
onClick=\"finalSubmit(%d, 'getSeqAlignment%d', 'getSeqGi', '%s%d', %d)\"><input \
type=\"hidden\" name=\"db\" value=\"\"><input type=\"hidden\" name=\"term\" \
value=\"\"><input type=\"hidden\" name=\"doptcmdl\" value=\"docsum\"><input \
type=\"hidden\" name=\"cmd\" value=\"search\"></form>",
                                     
                                     "<FORM  method=\"POST\" \
action=\"http://www.ncbi.nlm.nih.gov/Traces/trace.cgi\" \
name=\"%s%d\"><input type=button value=\"Get selected sequences\" \
onClick=\"finalSubmit(%d, 'getSeqAlignment%d', 'getSeqGi', '%s%d', %d)\"><input \
type=\"hidden\" name=\"val\" value=\"\"><input \
type=\"hidden\" name=\"cmd\" value=\"retrieve\"></form>"
};

static string k_GetSeqSelectForm = "<FORM><input \
type=\"button\" value=\"Select all\" onClick=\"handleCheckAll('select', \
'getSeqAlignment%d', 'getSeqGi')\"></form></td><td><FORM><input \
type=\"button\" value=\"Deselect all\" onClick=\"handleCheckAll('deselect', \
'getSeqAlignment%d', 'getSeqGi')\"></form>";


static string k_GetTreeViewForm =  "<FORM  method=\"post\" \
action=\"http://www.ncbi.nlm.nih.gov/blast/treeview/blast_tree_view.cgi?request=page&rid=%s&queryID=%s&distmode=on\" \
name=\"tree%s%d\" target=\"trv%s\"> \
<input type=button value=\"Distance tree of results\" onClick=\"extractCheckedSeq('getSeqAlignment%d', 'getSeqGi', 'tree%s%d')\"> \
<input type=\"hidden\" name=\"sequenceSet\" value=\"\"><input type=\"hidden\" name=\"screenWidth\" value=\"\"></form>";
#endif


static const int k_MaxDeflinesToShow = 8;
static const int k_MinDeflinesToShow = 3;


CDisplaySeqalign::CDisplaySeqalign(const CSeq_align_set& seqalign, 
                       CScope& scope,
                       list <CRef<CSeqLocInfo> >* mask_seqloc, 
                       list <FeatureInfo*>* external_feature,
                       const char* matrix_name /* = BLAST_DEFAULT_MATRIX */)
    : m_SeqalignSetRef(&seqalign),
      m_Seqloc(mask_seqloc),
      m_QueryFeature(external_feature),
      m_Scope(scope) 
{
    m_AlignOption = 0;
    m_SeqLocChar = eX;
    m_SeqLocColor = eBlack;
    m_LineLen = 60;
    m_IsDbNa = true;
    m_CanRetrieveSeq = false;
    m_DbName = NcbiEmptyString;
    m_NumAlignToShow = 1000000;
    m_AlignType = eNotSet;
    m_Rid = "0";
    m_CddRid = "0";
    m_EntrezTerm = NcbiEmptyString;
    m_QueryNumber = 0;
    m_BlastType = NcbiEmptyString;
    m_MidLineStyle = eBar;
    m_ConfigFile = NULL;
    m_Reg = NULL;
    m_DynamicFeature = NULL;    
    m_MasterGeneticCode = 1;
    m_SlaveGeneticCode = 1;
    m_AlignTemplates = NULL;
    m_Ctx = NULL;
    m_UseLinkoutDB = CLinkoutDB::UseLinkoutDB();

    CNcbiMatrix<int> mtx;
    CAlignFormatUtil::GetAsciiProteinMatrix(matrix_name 
                                       ? matrix_name 
                                       : BLAST_DEFAULT_MATRIX, mtx);
    _ASSERT(!mtx.GetData().empty());

    m_Matrix = new int*[mtx.GetRows()];
    for(size_t i = 0; i<mtx.GetRows(); ++i) {
        m_Matrix[i] = new int[mtx.GetCols()];
    }
    // copy data from matrix
    for(size_t i = 0; i<mtx.GetRows(); ++i) {
        for (size_t j = 0; j < mtx.GetCols(); j++) {
            m_Matrix[i][j] = mtx(i, j);
        }
    }
}


CDisplaySeqalign::~CDisplaySeqalign()
{
    for(int i = 0; i<k_NumAsciiChar; ++i) {
        delete [] m_Matrix[i];
    }
    delete [] m_Matrix;
    if (m_ConfigFile) {
        delete m_ConfigFile;
    } 
    if (m_Reg) {
        delete m_Reg;
    }
    
    if(m_DynamicFeature){
        delete m_DynamicFeature;
    }
}

//8.Display Identities,positives,frames etc
string CDisplaySeqalign::x_FormatIdentityInfo(string alignInfo, SAlnInfo* aln_vec_info)
{
    int aln_stop = (int)m_AV->GetAlnStop();       
    int master_strand  = m_AV->StrandSign(0);
    int slave_strand = m_AV->StrandSign(1);
    int master_frame = aln_vec_info->alnRowInfo->frame[0];
    int slave_frame = aln_vec_info->alnRowInfo->frame[1];
    bool aln_is_prot = (m_AlignType & eProt) != 0 ? true : false;


    string alignParams = alignInfo;//Some already filled in x_DisplayAlignInfo
    //out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
    //string alignParams = " Identities = " + NStr::IntToString(match);
    
    alignParams = CAlignFormatUtil::MapTemplate(alignParams, "aln_match",NStr::IntToString(aln_vec_info->match) + "/"+ NStr::IntToString(aln_stop+1));
    alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_ident",aln_vec_info->identity);
    
    if(aln_is_prot){
        //out<<", Positives = "<<(positive + match)<<"/"<<(aln_stop+1)
        //   <<" ("<<(((positive + match)*100)/(aln_stop+1))<<"%"<<")";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_pos",NStr::IntToString(aln_vec_info->positive + aln_vec_info->match) + "/" + NStr::IntToString(aln_stop+1));
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_pos_prc",NStr::IntToString(((aln_vec_info->positive + aln_vec_info->match)*100)/(aln_stop+1)));
    }
    else {//!!!!Check this!!!!
        //out<<" Strand="<<(master_strand==1 ? "Plus" : "Minus")
        //<<"/"<<(slave_strand==1? "Plus" : "Minus")<<"\n";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_strand",(master_strand==1 ? "Plus" : "Minus")+ (string)"/"+ (slave_strand==1? "Plus" : "Minus"));
    }
    //out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
   //    <<" ("<<((gap*100)/(aln_stop+1))<<"%"<<")"<<"\n";
    alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_gaps",NStr::IntToString(aln_vec_info->gap) + "/" + NStr::IntToString(aln_stop+1));
    alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_gaps_prc",NStr::IntToString((aln_vec_info->gap*100)/(aln_stop+1)));
    
    if(master_frame != 0 && slave_frame != 0) {
        //out <<" Frame = " << ((master_frame > 0) ? "+" : "") 
        //    << master_frame <<"/"<<((slave_frame > 0) ? "+" : "") 
        //    << slave_frame<<"\n";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame",((master_frame > 0) ? "+" : "") + master_frame + (string)"/"+((slave_frame > 0) ? "+" : ""));
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame_show","shown");        
    } else if (master_frame != 0){
        //out <<" Frame = " << ((master_frame > 0) ? "+" : "") 
        //    << master_frame << "\n";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame",(master_frame > 0) ? "+" : "");
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame_show","shown");        
    }  else if (slave_frame != 0){
        //out <<" Frame = " << ((slave_frame > 0) ? "+" : "") 
        //    << slave_frame <<"\n";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame",((slave_frame > 0) ? "+" : "") + slave_frame) ;
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame_show","shown");        
    }
    else {
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame","");        
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_frame_show","");        
    }     
    return alignParams;
}


///show blast identity, positive etc.
///@param out: output stream
///@param aln_stop: stop in aln coords
///@param identity: identity
///@param positive: positives
///@param match: match
///@param gap: gap
///@param master_strand: plus strand = 1 and minus strand = -1
///@param slave_strand:  plus strand = 1 and minus strand = -1
///@param master_frame: frame for master
///@param slave_frame: frame for slave
///@param aln_is_prot: is protein alignment?
///
static void s_DisplayIdentityInfo(CNcbiOstream& out, int aln_stop, 
                                   int identity, int positive, int match,
                                   int gap, int master_strand, 
                                   int slave_strand, int master_frame, 
                                   int slave_frame, bool aln_is_prot)
{
    out<<" Identities = "<<match<<"/"<<(aln_stop+1)<<" ("<<identity<<"%"<<")";
    if(aln_is_prot){
        out<<", Positives = "<<(positive + match)<<"/"<<(aln_stop+1)
           <<" ("<< CAlignFormatUtil::GetPercentMatch(positive + match, aln_stop+1) <<"%"<<")";
    }
    out<<", Gaps = "<<gap<<"/"<<(aln_stop+1)
       <<" ("<<((gap*100)/(aln_stop+1))<<"%"<<")"<<"\n";
    if (!aln_is_prot){ 
        out<<" Strand="<<(master_strand==1 ? "Plus" : "Minus")
           <<"/"<<(slave_strand==1? "Plus" : "Minus")<<"\n";
    }
    if(master_frame != 0 && slave_frame != 0) {
        out <<" Frame = " << ((master_frame > 0) ? "+" : "") 
            << master_frame <<"/"<<((slave_frame > 0) ? "+" : "") 
            << slave_frame<<"\n";
    } else if (master_frame != 0){
        out <<" Frame = " << ((master_frame > 0) ? "+" : "") 
            << master_frame << "\n";
    }  else if (slave_frame != 0){
        out <<" Frame = " << ((slave_frame > 0) ? "+" : "") 
            << slave_frame <<"\n";
    } 
    out<<"\n";
    
}

///wrap line
///@param out: output stream
///@param str: string to wrap
///
static void s_WrapOutputLine(CNcbiOstream& out, const string& str)
{
    const int line_len = 60;
    bool do_wrap = false;
    int length = (int) str.size();
    if (length > line_len) {
        for (int i = 0; i < length; i ++){
            if(i > 0 && i % line_len == 0){
                do_wrap = true;
            }   
            out << str[i];
            if(do_wrap && isspace((unsigned char) str[i])){
                out << "\n";  
                do_wrap = false;
            }
        }
    } else {
        out << str;
    }
}

///To add color to bases other than identityChar
///@param seq: sequence
///@param identity_char: identity character
///@param out: output stream
///
static void s_ColorDifferentBases(string& seq, char identity_char,
                                  CNcbiOstream& out){
    string base_color = k_ColorRed;
    bool tagOpened = false;
    for(int i = 0; i < (int)seq.size(); i ++){
        if(seq[i] != identity_char){
            if(!tagOpened){
                out << "<font color=\""+base_color+"\"><b>";
                tagOpened =  true;
            }
            
        } else {
            if(tagOpened){
                out << "</b></font>";
                tagOpened = false;
            }
        }
        out << seq[i];
        if(tagOpened && i == (int)seq.size() - 1){
            out << "</b></font>";
            tagOpened = false;
        }
    } 
}

///return the frame for a given strand
///Note that start is zero bases.  It returns frame +/-(1-3). 0 indicates error
///@param start: sequence start position
///@param strand: strand
///@param id: the seqid
///@param scope: the scope
///@return: the frame
///
static int s_GetFrame (int start, ENa_strand strand, const CSeq_id& id, 
                       CScope& sp) 
{
    int frame = 0;
    if (strand == eNa_strand_plus) {
        frame = (start % 3) + 1;
    } else if (strand == eNa_strand_minus) {
        frame = -(((int)sp.GetBioseqHandle(id).GetBioseqLength() - start - 1)
                  % 3 + 1);
        
    }
    return frame;
}

///reture the frame for master seq in stdseg
///@param ss: the input stdseg
///@param scope: the scope
///@return: the frame
///
static int s_GetStdsegMasterFrame(const CStd_seg& ss, CScope& scope)
{
    const CRef<CSeq_loc> slc = ss.GetLoc().front();
    ENa_strand strand = GetStrand(*slc);
    int frame = s_GetFrame(strand ==  eNa_strand_plus ?
                           GetStart(*slc, &scope) : GetStop(*slc, &scope),
                           strand ==  eNa_strand_plus ?
                           eNa_strand_plus : eNa_strand_minus,
                           *(ss.GetIds().front()), scope);
    return frame;
}

///return the get sequence table for html display
///@param form_name: form name
///@parm db_is_na: is the db of nucleotide type?
///@query_number: the query number
///@return: the from string
///
static string s_GetSeqForm(char* form_name, bool db_is_na, int query_number,
                           int db_type, const string& dbName,const char *rid, const char *queryID,bool showTreeButtons)
{
    string temp_buf = NcbiEmptyString;
    AutoPtr<char, ArrayDeleter<char> > buf(new char[dbName.size() + 4096]);
    if(form_name){             
        string localClientButtons = "";
        if(showTreeButtons) {
            string l_GetTreeViewForm  = CAlignFormatUtil::GetURLFromRegistry( "TREEVIEW_FRM");
            localClientButtons = "<td>" + l_GetTreeViewForm + "</td>";
        }
        string l_GetSeqSubmitForm  = CAlignFormatUtil::GetURLFromRegistry( "GETSEQ_SUB_FRM", db_type); 
        string l_GetSeqSelectForm = CAlignFormatUtil::GetURLFromRegistry( "GETSEQ_SEL_FRM");
        
        string template_str = "<table border=\"0\"><tr><td>" +
                            l_GetSeqSubmitForm + //k_GetSeqSubmitForm[db_type] +
                            "</td><td>" +
                             l_GetSeqSelectForm + //k_GetSeqSelectForm +
                            "</td>" +
                            localClientButtons +
                            "</tr></table>";

        if(showTreeButtons) {
            sprintf(buf.get(), template_str.c_str(), form_name, query_number,
                db_is_na?1:0, query_number, form_name, query_number, db_type, 
                query_number,query_number,             
                    rid,queryID,form_name,query_number,rid,query_number,form_name,query_number);                  
        
        }
        else {
             sprintf(buf.get(), template_str.c_str(), form_name, query_number,
                db_is_na?1:0, query_number, form_name, query_number, db_type, 
                query_number,query_number);              
        }

    }
    temp_buf = buf.get();
    return temp_buf;
}

///Gets Query Seq ID from Seq Align
///@param actual_aln_list: align set for one query
///@return: string query ID
static string s_GetQueryIDFromSeqAlign(const CSeq_align_set& actual_aln_list) 
{
    CRef<CSeq_align> first_aln = actual_aln_list.Get().front();
    const CSeq_id& query_SeqID = first_aln->GetSeq_id(0);
    string queryID;
    query_SeqID.GetLabel(&queryID);    
    return queryID;
}

///return id type specified or null ref
///@param ids: the input ids
///@param choice: id of choice
///@return: the id with specified type
///
static CRef<CSeq_id> s_GetSeqIdByType(const list<CRef<CSeq_id> >& ids, 
                                      CSeq_id::E_Choice choice)
{
    CRef<CSeq_id> cid;
    
    for (CBioseq::TId::const_iterator iter = ids.begin(); iter != ids.end(); 
         iter ++){
        if ((*iter)->Which() == choice){
            cid = *iter;
            break;
        }
    }
    
    return cid;
}

///return gi from id list
///@param ids: the input ids
///@return: the gi if found
///
static int s_GetGiForSeqIdList (const list<CRef<CSeq_id> >& ids)
{
    int gi = 0;
    CRef<CSeq_id> id = s_GetSeqIdByType(ids, CSeq_id::e_Gi);
    if (!(id.Empty())){
        return id->GetGi();
    }
    return gi;
}


///return concatenated exon sequence
///@param feat: the feature containing this cds
///@param feat_strand: the feature strand
///@param range: the range list of seqloc
///@param total_coding_len: the total exon length excluding intron
///@param raw_cdr_product: the raw protein sequence
///@return: the concatenated exon sequences with amino acid aligned to 
///to the second base of a codon
///
static string s_GetConcatenatedExon(CFeat_CI& feat,  
                                    ENa_strand feat_strand, 
                                    list<CRange<TSeqPos> >& range,
                                    TSeqPos total_coding_len,
                                    string& raw_cdr_product, TSeqPos frame_adj)
{

    string concat_exon(total_coding_len, ' ');
    TSeqPos frame = 1;
    const CCdregion& cdr = feat->GetData().GetCdregion();
    if(cdr.IsSetFrame()){
        frame = cdr.GetFrame();
    }
    TSeqPos num_coding_base;
    int num_base;
    TSeqPos coding_start_base;
    if(feat_strand == eNa_strand_minus){
        coding_start_base = total_coding_len - 1 - (frame -1) - frame_adj;
        num_base = total_coding_len - 1;
        num_coding_base = 0;
        
    } else {
        coding_start_base = 0; 
        coding_start_base += frame - 1 + frame_adj;
        num_base = 0;
        num_coding_base = 0;
    }
    
    ITERATE(list<CRange<TSeqPos> >, iter, range){
        //note that feature on minus strand needs to be
        //filled backward.
        if(feat_strand != eNa_strand_minus){
            for(TSeqPos i = 0; i < iter->GetLength(); i ++){
                if((TSeqPos)num_base >= coding_start_base){
                    num_coding_base ++;
                    if(num_coding_base % 3 == 2){
                        //a.a to the 2nd base
                        if(num_coding_base / 3 < raw_cdr_product.size()){
                            //make sure the coding region is no
                            //more than the protein seq as there
                            //could errors in ncbi record
                            concat_exon[num_base] 
                                = raw_cdr_product[num_coding_base / 3];
                        }                           
                    }
                }
                num_base ++;
            }    
        } else {
            
            for(TSeqPos i = 0; i < iter->GetLength() &&
                    num_base >= 0; i ++){
                if((TSeqPos)num_base <= coding_start_base){
                    num_coding_base ++;
                    if(num_coding_base % 3 == 2){
                        //a.a to the 2nd base
                        if(num_coding_base / 3 < 
                           raw_cdr_product.size() &&
                           coding_start_base >= num_coding_base){
                            //make sure the coding region is no
                            //more than the protein seq as there
                            //could errors in ncbi record
                            concat_exon[num_base] 
                                = raw_cdr_product[num_coding_base / 3];
                        }                           
                    }
                }
                num_base --;
            }    
        }
    }
    return concat_exon;
}

///map slave feature info to master seq
///@param master_feat_range: master feature seqloc to be filled
///@param feat: the feature in concern
///@param slave_feat_range: feature info for slave
///@param av: the alignment vector for master-slave seqalign
///@param row: the row
///@param frame_adj: frame adjustment
///

static void s_MapSlaveFeatureToMaster(list<CRange<TSeqPos> >& master_feat_range,
                                      ENa_strand& master_feat_strand, CFeat_CI& feat,
                                      list<CSeq_loc_CI::TRange>& slave_feat_range,
                                      ENa_strand slave_feat_strand,
                                      CAlnVec* av, 
                                      int row, TSeqPos frame_adj)
{
    TSeqPos trans_frame = 1;
    const CCdregion& cdr = feat->GetData().GetCdregion();
    if(cdr.IsSetFrame()){
        trans_frame = cdr.GetFrame();
    }
    trans_frame += frame_adj;

    TSeqPos prev_exon_len = 0;
    bool is_first_in_range = true;

    if ((av->IsPositiveStrand(1) && slave_feat_strand == eNa_strand_plus) || 
        (av->IsNegativeStrand(1) && slave_feat_strand == eNa_strand_minus)) {
        master_feat_strand = eNa_strand_plus;
    } else {
        master_feat_strand = eNa_strand_minus;
    }
    
    list<CSeq_loc_CI::TRange> acutal_slave_feat_range = slave_feat_range;

    ITERATE(list<CSeq_loc_CI::TRange>, iter_temp,
            acutal_slave_feat_range){
        CRange<TSeqPos> actual_feat_seq_range = av->GetSeqRange(row).
            IntersectionWith(*iter_temp);          
        if(!actual_feat_seq_range.Empty()){
            TSeqPos slave_aln_from = 0, slave_aln_to = 0;
            TSeqPos frame_offset = 0;
            int curr_exon_leading_len = 0;
            //adjust frame 
            if (is_first_in_range) {  
                if (slave_feat_strand == eNa_strand_plus) {
                    curr_exon_leading_len 
                        = actual_feat_seq_range.GetFrom() - iter_temp->GetFrom();
                    
                } else {
                    curr_exon_leading_len 
                        = iter_temp->GetTo() - actual_feat_seq_range.GetTo();
                }
                is_first_in_range = false;
                frame_offset = (3 - (prev_exon_len + curr_exon_leading_len)%3
                                + (trans_frame - 1)) % 3;
            }
            
            if (av->IsPositiveStrand(1) && 
                slave_feat_strand == eNa_strand_plus) {
                slave_aln_from 
                    = av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetFrom() + 
                                              frame_offset, CAlnMap::eRight );   
                
                slave_aln_to =
                    av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetTo(),
                                            CAlnMap::eLeft);
            } else if (av->IsNegativeStrand(1) && 
                       slave_feat_strand == eNa_strand_plus) {
                
                slave_aln_from 
                    = av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetTo(),
                                              CAlnMap::eRight);   
                
                slave_aln_to =
                    av->GetAlnPosFromSeqPos(row,
                                            actual_feat_seq_range.GetFrom() +
                                            frame_offset, CAlnMap::eLeft);
            }  else if (av->IsPositiveStrand(1) && 
                        slave_feat_strand == eNa_strand_minus) {
                slave_aln_from 
                    = av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetFrom(),
                                              CAlnMap::eRight);   
                
                slave_aln_to =
                    av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetTo() -
                                            frame_offset, CAlnMap::eLeft);

            } else if (av->IsNegativeStrand(1) && 
                       slave_feat_strand == eNa_strand_minus){
                slave_aln_from 
                    = av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetTo() - 
                                              frame_offset, CAlnMap::eRight );   
                
                slave_aln_to =
                    av->GetAlnPosFromSeqPos(row, actual_feat_seq_range.GetFrom(),
                                            CAlnMap::eLeft);
            }
            
            TSeqPos master_from = 
                av->GetSeqPosFromAlnPos(0, slave_aln_from, CAlnMap::eRight);
            
            TSeqPos master_to = 
                av->GetSeqPosFromAlnPos(0, slave_aln_to, CAlnMap::eLeft);
            
            CRange<TSeqPos> master_range(master_from, master_to);
            master_feat_range.push_back(master_range); 
            
        }
        prev_exon_len += iter_temp->GetLength();
    }
}



///return cds coded sequence and fill the id if found
///@param genetic_code: the genetic code
///@param feat: the feature containing this cds
///@param scope: scope to fetch sequence
///@param range: the range list of seqloc
///@param handle: the bioseq handle
///@param feat_strand: the feature strand
///@param feat_id: the feature id to be filled
///@param frame_adj: frame adjustment
///@param mix_loc: is this seqloc mixed with other seqid?
///@return: the encoded protein sequence
///
static string s_GetCdsSequence(int genetic_code, CFeat_CI& feat, 
                               CScope& scope, list<CRange<TSeqPos> >& range,
                               const CBioseq_Handle& handle, 
                               ENa_strand feat_strand, string& feat_id,
                               TSeqPos frame_adj, bool mix_loc)
{
    string raw_cdr_product = NcbiEmptyString;
    if(feat->IsSetProduct() && feat->GetProduct().IsWhole() && !mix_loc){
        //show actual aa  if there is a cds product
          
        const CSeq_id& productId = 
            feat->GetProduct().GetWhole();
        const CBioseq_Handle& productHandle 
            = scope.GetBioseqHandle(productId );
        feat_id = "CDS:" + 
            GetTitle(productHandle).substr(0, k_FeatureIdLen);
        productHandle.
            GetSeqVector(CBioseq_Handle::eCoding_Iupac).
            GetSeqData(0, productHandle.
            GetBioseqLength(), raw_cdr_product);
    } else { 
        CSeq_loc isolated_loc;
        ITERATE(list<CRange<TSeqPos> >, iter, range){
            TSeqPos from = iter->GetFrom();
            TSeqPos to = iter->GetTo();
            if(feat_strand == eNa_strand_plus){
                isolated_loc.
                    Add(*(handle.GetRangeSeq_loc(from + frame_adj,
                                                 to,         
                                                 feat_strand)));
            } else {
                isolated_loc.
                    Add(*(handle.GetRangeSeq_loc(from,
                                                 to - frame_adj,   
                                                 feat_strand)));
            }
        }
        CGenetic_code gc;
        CRef<CGenetic_code::C_E> ce(new CGenetic_code::C_E);
        ce->Select(CGenetic_code::C_E::e_Id);
        ce->SetId(genetic_code);
        gc.Set().push_back(ce);
        CSeqTranslator::Translate(isolated_loc, handle.GetScope(),
                                  raw_cdr_product, &gc);
      
    }
    return raw_cdr_product;
}

///fill the cds start positions (1 based)
///@param line: the input cds line
///@param concat_exon: exon only string
///@param length_per_line: alignment length per line
///@param feat_aln_start_totalexon: feature aln pos in concat_exon
///@param strand: the alignment strand
///@param start: start list to be filled
///
static void s_FillCdsStartPosition(string& line, string& concat_exon,
                                   size_t length_per_line,
                                   TSeqPos feat_aln_start_totalexon,
                                   ENa_strand seq_strand, 
                                   ENa_strand feat_strand,
                                   list<TSeqPos>& start)
{
    size_t actual_line_len = 0;
    size_t aln_len = line.size();
    TSeqPos previous_num_letter = 0;
    
    //the number of amino acids preceeding this exon start position
    for (size_t i = 0; i <= feat_aln_start_totalexon; i ++){
        if(feat_strand == eNa_strand_minus){
            //remember the amino acid in this case goes backward
            //therefore we count backward too
            
            int pos = concat_exon.size() -1 - i;
            if(pos >= 0 && isalpha((unsigned char) concat_exon[pos])){
                previous_num_letter ++;
            }
            
        } else {
            if(isalpha((unsigned char) concat_exon[i])){
                previous_num_letter ++;
            }
        }
    }
    
    
    TSeqPos prev_num = 0;
    //go through the entire feature line and get the amino acid position 
    //for each line
    for(size_t i = 0; i < aln_len; i += actual_line_len){
        //handle the last row which may be shorter
        if(aln_len - i< length_per_line) {
            actual_line_len = aln_len - i;
        } else {
            actual_line_len = length_per_line;
        }
        //the number of amino acids on this row
        TSeqPos cur_num = 0;
        bool has_intron = false;
        
        //go through each character on a row
        for(size_t j = i; j < actual_line_len + i; j ++){
            //don't count gap
            if(isalpha((unsigned char) line[j])){
                cur_num ++;
            } else if(line[j] == k_IntronChar){
                has_intron = true;
            }
        }
            
        if(cur_num > 0){
            if(seq_strand == eNa_strand_plus){
                if(feat_strand == eNa_strand_minus) {
                    start.push_back(previous_num_letter - prev_num); 
                } else {
                    start.push_back(previous_num_letter + prev_num);  
                }
            } else {
                if(feat_strand == eNa_strand_minus) {
                    start.push_back(previous_num_letter + prev_num);  
                } else {
                    start.push_back(previous_num_letter - prev_num);  
                } 
            }
        } else if (has_intron) {
            start.push_back(0);  //sentinal for no show
        }
        prev_num += cur_num;
    }
}

///make a new copy of master seq with feature info and return the scope
///that contains this sequence
///@param feat_range: the feature seqlocs
///@param feat_seq_strand: the stand info
///@param handle: the seq handle for the original master seq
///@return: the scope containing the new master seq
///
static CRef<CScope> s_MakeNewMasterSeq(list<list<CRange<TSeqPos> > >& feat_range,
                                       list<ENa_strand>& feat_seq_strand,
                                       const CBioseq_Handle& handle) 
{
    CRef<CObjectManager> obj;
    obj = CObjectManager::GetInstance();
    CGBDataLoader::RegisterInObjectManager(*obj);       
    CRef<CScope> scope (new CScope(*obj));
    scope->AddDefaults();
    CRef<CBioseq> cbsp(new CBioseq());
    cbsp->Assign(*(handle.GetCompleteBioseq()));

    CBioseq::TAnnot& anot_list = cbsp->SetAnnot();
    CRef<CSeq_annot> anot(new CSeq_annot);
    CRef<CSeq_annot::TData> data(new CSeq_annot::TData);
    data->Select(CSeq_annot::TData::e_Ftable);
    anot->SetData(*data);
    CSeq_annot::TData::TFtable& ftable = anot->SetData().SetFtable();
    int counter = 0;
    ITERATE(list<list<CRange<TSeqPos> > >, iter, feat_range) {
        counter ++;
        CRef<CSeq_feat> seq_feat(new CSeq_feat);
        CRef<CSeqFeatData> feat_data(new CSeqFeatData);
        feat_data->Select(CSeq_feat::TData::e_Cdregion);
        seq_feat->SetData(*feat_data);
        seq_feat->SetComment("Putative " + NStr::IntToString(counter));
        CRef<CSeq_loc> seq_loc (new CSeq_loc);
       
        ITERATE(list<CRange<TSeqPos> >, iter2, *iter) {
            seq_loc->Add(*(handle.GetRangeSeq_loc(iter2->GetFrom(),
                                                  iter2->GetTo(),
                                                  feat_seq_strand.front())));
        }
        seq_feat->SetLocation(*seq_loc);
        ftable.push_back(seq_feat);
        feat_seq_strand.pop_front();
    }
    anot_list.push_back(anot);
    CRef<CSeq_entry> entry(new CSeq_entry());
    entry->SetSeq(*cbsp);
    scope->AddTopLevelSeqEntry(*entry);
  
    return scope;
}

//output feature lines
//@param reference_feat_line: the master feature line to be compared 
//for coloring
//@param feat_line: the slave feature line
//@param color_feat_mismatch: color or not
//@param start: the alignment pos
//@param len: the length per line
//@param out: stream for output
//
static void s_OutputFeature(string& reference_feat_line, 
                            string& feat_line,
                            bool color_feat_mismatch,
                            int start,
                            int len,  
                            CNcbiOstream& out)
{
    if((int)feat_line.size() > start){
        string actual_feat = feat_line.substr(start, len);
        string actual_reference_feat = NcbiEmptyString;
        if(reference_feat_line != NcbiEmptyString){
            actual_reference_feat = reference_feat_line.substr(start, len);
        }
        if(color_feat_mismatch 
           && actual_reference_feat != NcbiEmptyString &&
           !NStr::IsBlank(actual_reference_feat)){
            string base_color = k_ColorPink;
            bool tagOpened = false;
            for(int i = 0; i < (int)actual_feat.size() &&
                    i < (int)actual_reference_feat.size(); i ++){
                if (actual_feat[i] != actual_reference_feat[i]) {
                    if(actual_feat[i] != ' ' &&
                       actual_feat[i] != k_IntronChar &&
                       actual_reference_feat[i] != k_IntronChar) {
                        if(!tagOpened){
                            out << "<font color=\""+base_color+"\"><b>";
                            tagOpened =  true;
                        }
                        
                    }
                } else {
                    if (actual_feat[i] != ' '){ //no close if space to 
                        //minimizing the open and close of tags
                        if(tagOpened){
                            out << "</b></font>";
                            tagOpened = false;
                        }
                    }
                }
                out << actual_feat[i];
                //close tag at the end of line
                if(tagOpened && i == (int)actual_feat.size() - 1){
                    out << "</b></font>";
                    tagOpened = false;
                }
            }
        } else {
            out << actual_feat;
        }
    }
    
}


void CDisplaySeqalign::x_PrintFeatures(TSAlnFeatureInfoList& feature,
                                       int row, 
                                       CAlnMap::TSignedRange alignment_range,
                                       int aln_start,
                                       int line_length, 
                                       int id_length,
                                       int start_length,
                                       int max_feature_num, 
                                       string& master_feat_str,
                                       CNcbiOstream& out)
{
    NON_CONST_ITERATE(TSAlnFeatureInfoList, iter, feature) {
        //check blank string for cases where CDS is in range 
        //but since it must align with the 2nd codon and is 
        //actually not in range
        if (alignment_range.IntersectingWith((*iter)->aln_range) && 
            !(NStr::IsBlank((*iter)->feature_string.
                            substr(aln_start, line_length)) &&
              m_AlignOption & eShowCdsFeature)){  
            if((m_AlignOption&eHtml)&&(m_AlignOption&eMergeAlign)
               && (m_AlignOption&eSequenceRetrieval && m_CanRetrieveSeq)){
                char checkboxBuf[200];
                sprintf(checkboxBuf,  k_UncheckabeCheckbox.c_str(),
                        m_QueryNumber);
                out << checkboxBuf;
            }
            out<<(*iter)->feature->feature_id;
            if((*iter)->feature_start.empty()){
                CAlignFormatUtil::
                    AddSpace(out, id_length + k_IdStartMargin
                             +start_length + k_StartSequenceMargin
                             -(*iter)->feature->feature_id.size());
            } else {
                int feat_start = (*iter)->feature_start.front();
                if(feat_start > 0){
                    CAlignFormatUtil::
                        AddSpace(out, id_length + k_IdStartMargin
                                 -(*iter)->feature->feature_id.size());
                    out << feat_start;
                    CAlignFormatUtil::
                        AddSpace(out, start_length -
                                 NStr::IntToString(feat_start).size() +
                                 k_StartSequenceMargin);
                } else { //no show start
                    CAlignFormatUtil::
                        AddSpace(out, id_length + k_IdStartMargin
                                 +start_length + k_StartSequenceMargin
                                 -(*iter)->feature->feature_id.size());
                }
                
                (*iter)->feature_start.pop_front();
            }
            bool color_cds_mismatch = false; 
            if(max_feature_num == 1 && (m_AlignOption & eHtml) && 
               (m_AlignOption & eShowCdsFeature) && row > 0){
                //only for slaves, only for cds feature
                //only color mismach if only one cds exists
                color_cds_mismatch = true;
            }
            s_OutputFeature(master_feat_str, 
                            (*iter)->feature_string,
                            color_cds_mismatch, aln_start,
                            line_length, out);
            if(row == 0){//set master feature as reference
                master_feat_str = (*iter)->feature_string;
            }
            out<<"\n";
        }
    }
    
}



string CDisplaySeqalign::x_GetUrl(CAlignFormatUtil::SSeqURLInfo *seqUrlInfo, const list<CRef<CSeq_id> >& ids,bool useTemplates) const
{
    string urlLink = NcbiEmptyString;    
 
    seqUrlInfo->gi = (seqUrlInfo->gi == 0) ? s_GetGiForSeqIdList(ids):seqUrlInfo->gi;
    seqUrlInfo->user_url= m_Reg->Get(m_BlastType, "TOOL_URL");
    
    urlLink = CAlignFormatUtil::GetIDUrl(seqUrlInfo,&ids,useTemplates);
    return urlLink;
}



void
CDisplaySeqalign::SetSubjectMasks(const TSeqLocInfoVector& masks)
{
    ITERATE(TSeqLocInfoVector, sequence_masks, masks) {
        const CSeq_id& id = sequence_masks->front()->GetSeqId();
        m_SubjectMasks[id] = *sequence_masks;
    }
}

CDisplaySeqalign::SAlnRowInfo *CDisplaySeqalign::x_PrepareRowData(void)
{
    size_t maxIdLen=0, maxStartLen=0;
    //, startLen=0, actualLineLen=0;
    //size_t aln_stop=m_AV->GetAlnStop();
    const int rowNum=m_AV->GetNumRows();   
    if(m_AlignOption & eMasterAnchored){
        m_AV->SetAnchor(0);
    }
    m_AV->SetGapChar('-');

    if (m_AlignOption & eShowEndGaps) {
        m_AV->SetEndChar('-');
    }
    else {
        m_AV->SetEndChar(' ');
    }    
    vector<string> sequence(rowNum);
    vector<CAlnMap::TSeqPosList> seqStarts(rowNum);
    vector<CAlnMap::TSeqPosList> seqStops(rowNum);
    vector<CAlnMap::TSeqPosList> insertStart(rowNum);
    vector<CAlnMap::TSeqPosList> insertAlnStart(rowNum);
    vector<CAlnMap::TSeqPosList> insertLength(rowNum);
    vector<string> seqidArray(rowNum);
    string middleLine;
    vector<CAlnMap::TSignedRange> rowRng(rowNum);
    vector<int> frame(rowNum);
    vector<int> taxid(rowNum);
    int max_feature_num = 0;
    
    
    //Add external query feature info such as phi blast pattern
    vector<TSAlnFeatureInfoList> bioseqFeature;
    x_GetQueryFeatureList(rowNum, (int)m_AV->GetAlnStop(), bioseqFeature);
    _ASSERT((int)bioseqFeature.size() == rowNum);
    // Mask locations for queries (first elem) and subjects (all other rows)
    vector<TSAlnSeqlocInfoList> masked_regions(rowNum);
    x_FillLocList(masked_regions[0], m_Seqloc);

    for (int row = 1; row < rowNum; row++) {
        const CSeq_id& id = m_AV->GetSeqId(row);
        x_FillLocList(masked_regions[row], &m_SubjectMasks[id]);
    }
    
    //prepare data for each row 
    list<list<CRange<TSeqPos> > > feat_seq_range;
    list<ENa_strand> feat_seq_strand;
    for (int row=0; row<rowNum; row++) {
        string type_temp = m_BlastType;
        type_temp = NStr::TruncateSpaces(NStr::ToLower(type_temp));
        if(type_temp == "mapview" || type_temp == "mapview_prev" || 
           type_temp == "gsfasta" || type_temp == "gsfasta_prev"){
            taxid[row] = CAlignFormatUtil::GetTaxidForSeqid(m_AV->GetSeqId(row),
                                                            m_Scope);
        } else {
            taxid[row] = 0;
        }
        rowRng[row] = m_AV->GetSeqAlnRange(row);
        frame[row] = (m_AV->GetWidth(row) == 3 ? 
                      s_GetFrame(m_AV->IsPositiveStrand(row) ? 
                                 m_AV->GetSeqStart(row) : 
                                 m_AV->GetSeqStop(row), 
                                 m_AV->IsPositiveStrand(row) ? 
                                 eNa_strand_plus : eNa_strand_minus, 
                                 m_AV->GetSeqId(row), m_Scope) : 0);        
        //make sequence
        m_AV->GetWholeAlnSeqString(row, sequence[row], &insertAlnStart[row],
                                   &insertStart[row], &insertLength[row],
                                   (int)m_LineLen, &seqStarts[row], &seqStops[row]);
        //make feature. Only for pairwise and untranslated for subject nuc seq
        if(!(m_AlignOption & eMasterAnchored) &&
           !(m_AlignOption & eMergeAlign) && m_AV->GetWidth(row) != 3 &&
           !(m_AlignType & eProt)){
            if(m_AlignOption & eShowCdsFeature){
                int master_gi = FindGi(m_AV->GetBioseqHandle(0).
                                       GetBioseqCore()->GetId());
                x_GetFeatureInfo(bioseqFeature[row], *m_featScope, 
                                 CSeqFeatData::e_Cdregion, row, sequence[row],
                                 feat_seq_range, feat_seq_strand,
                                 row == 1 && !(master_gi > 0) ? true : false);
                
                if(!(feat_seq_range.empty()) && row == 1) {
                    //make a new copy of master bioseq and add the feature from
                    //slave to make putative cds feature 
                    CRef<CScope> master_scope_with_feat = 
                        s_MakeNewMasterSeq(feat_seq_range, feat_seq_strand,
                                           m_AV->GetBioseqHandle(0));
                    //make feature string for master bioseq
                    list<list<CRange<TSeqPos> > > temp_holder;
                    x_GetFeatureInfo(bioseqFeature[0], *master_scope_with_feat, 
                                     CSeqFeatData::e_Cdregion, 0, sequence[0],
                                     temp_holder, feat_seq_strand, false);
                }
            }
            if(m_AlignOption & eShowGeneFeature){
                x_GetFeatureInfo(bioseqFeature[row], *m_featScope,
                                 CSeqFeatData::e_Gene, row, sequence[row],
                                 feat_seq_range, feat_seq_strand, false);
            }
        }
        //make id
        x_FillSeqid(seqidArray[row], row);
        maxIdLen=max<size_t>(seqidArray[row].size(), maxIdLen);
        size_t maxCood=max<size_t>(m_AV->GetSeqStart(row), m_AV->GetSeqStop(row));
        maxStartLen = max<size_t>(NStr::IntToString(maxCood).size(), maxStartLen);
    }
    for(int i = 0; i < rowNum; i ++){//adjust max id length for feature id 
        int num_feature = 0;
        ITERATE(TSAlnFeatureInfoList, iter, bioseqFeature[i]) {
            maxIdLen=max<size_t>((*iter)->feature->feature_id.size(), maxIdLen);
            num_feature ++;
            if(num_feature > max_feature_num){
                max_feature_num = num_feature;
            }
        }
    }  //end of preparing row data
    SAlnRowInfo *alnRoInfo = new SAlnRowInfo();
    alnRoInfo->sequence = sequence;
    alnRoInfo->seqStarts = seqStarts;
    alnRoInfo->seqStops = seqStops;
    alnRoInfo->insertStart = insertStart;
    alnRoInfo->insertAlnStart = insertAlnStart;
    alnRoInfo->insertLength = insertLength;
    alnRoInfo->seqidArray = seqidArray;
    alnRoInfo->middleLine = middleLine;
    alnRoInfo->rowRng = rowRng;
    alnRoInfo->frame = frame;
    alnRoInfo->taxid = taxid;
    alnRoInfo->bioseqFeature = bioseqFeature;
    alnRoInfo->masked_regions = masked_regions;
    alnRoInfo->seqidArray = seqidArray;
    alnRoInfo->maxIdLen = maxIdLen;
    alnRoInfo->maxStartLen = maxStartLen;
    alnRoInfo->max_feature_num = max_feature_num;    
    alnRoInfo->colorMismatch = false;
    alnRoInfo->rowNum = rowNum;
    return alnRoInfo;
}
//uses m_AV    m_LineLen m_AlignOption m_QueryNumber
string CDisplaySeqalign::x_DisplayRowData(SAlnRowInfo *alnRoInfo)
{
    size_t startLen=0, actualLineLen=0;
    string master_feat_str = NcbiEmptyString;
    size_t aln_stop=m_AV->GetAlnStop();

    int rowNum = alnRoInfo->rowNum;
    vector<int> prev_stop(rowNum);
    CNcbiOstrstream out;
    //output rows    
    for(int j=0; j<=(int)aln_stop; j+=(int)m_LineLen){
        //output according to aln coordinates
        if(aln_stop-j+1<m_LineLen) {
            actualLineLen=aln_stop-j+1;
        } else {
            actualLineLen=m_LineLen;
        }
        CAlnMap::TSignedRange curRange(j, j+(int)actualLineLen-1);
        //here is each row
        for (int row=0; row<rowNum; row++) {
            bool hasSequence = true;   
            if (!(m_AlignOption & eShowGapOnlyLines)) {
                hasSequence = curRange.IntersectingWith(alnRoInfo->rowRng[row]);
            }
            //only output rows that have sequence
            if (hasSequence){
                int start = alnRoInfo->seqStarts[row].front() + 1;  //+1 for 1 based
                int end = alnRoInfo->seqStops[row].front() + 1;
                list<string> inserts;
                string insertPosString;  //the one with "\" to indicate insert
                if(m_AlignOption & eMasterAnchored){
                    TSInsertInformationList insertList;
                    x_GetInserts(insertList, alnRoInfo->insertAlnStart[row], 
                                 alnRoInfo->insertStart[row], alnRoInfo->insertLength[row],  
                                 j + (int)m_LineLen);
                    x_FillInserts(row, curRange, j, inserts, insertPosString, 
                                  insertList);
                }
                //feature for query
                if(row == 0){          
                    x_PrintFeatures(alnRoInfo->bioseqFeature[row], row, curRange,
                                    j,(int)actualLineLen, alnRoInfo->maxIdLen,
                                    alnRoInfo->maxStartLen, alnRoInfo->max_feature_num, 
                                    master_feat_str, out); 
                }

                string urlLink = NcbiEmptyString;
                //setup url link for seqid
                if(m_AlignOption & eHtml){
                    int gi = 0;
                    if(m_AV->GetSeqId(row).Which() == CSeq_id::e_Gi){
                        gi = m_AV->GetSeqId(row).GetGi();
                    }
                    if(!(gi > 0)){
                        gi = s_GetGiForSeqIdList(m_AV->GetBioseqHandle(row).
                                                 GetBioseqCore()->GetId());
                    }
                    if((row == 0 && (m_AlignOption & eHyperLinkMasterSeqid)) ||
                       (row > 0 && (m_AlignOption & eHyperLinkSlaveSeqid))){
                        
                        if(gi > 0){
                            out<<"<a name="<<gi<<"></a>";
                        } else {
                            out<<"<a name="<<alnRoInfo->seqidArray[row]<<"></a>";
                        }
                    }
                    //get sequence checkbox
                    if((m_AlignOption & eMergeAlign) && 
                        (m_AlignOption & eSequenceRetrieval) && m_CanRetrieveSeq){
                        char checkboxBuf[512];
                        if (row == 0) {
                            sprintf(checkboxBuf, k_UncheckabeCheckbox.c_str(),
                                    m_QueryNumber); 
                        } else {
                            sprintf(checkboxBuf, k_Checkbox.c_str(), gi > 0 ?
                                    NStr::IntToString(gi).c_str() :
                                    alnRoInfo->seqidArray[row].c_str(), m_QueryNumber);
                        }
                        out << checkboxBuf;        
                    }
                    else if(m_AlignOption & eShowCheckBox) {                        
                        const CRef<CSeq_id> seqID = FindBestChoice(m_AV->GetBioseqHandle(row).GetBioseqCore()->GetId(), CSeq_id::WorstRank);
                        string id_str = CAlignFormatUtil::GetLabel(seqID);
                        if(seqID->IsLocal()) {
                            id_str = "lcl|" + id_str;            
                        }        
                        char checkboxBuf[512];                        
                        sprintf(checkboxBuf, k_CheckboxEx.c_str(), id_str.c_str());
                        out << checkboxBuf;        
                    }                    
                    if((row == 0 && (m_AlignOption & eHyperLinkMasterSeqid)) ||
                       (row > 0 && (m_AlignOption & eHyperLinkSlaveSeqid))){
                        string idString = m_AV->GetSeqId(row).GetSeqIdString();												
                        CRange<TSeqPos> range = (m_AlnLinksParams.count(idString) > 0 && m_AlnLinksParams[idString].subjRange) ?
								CRange<TSeqPos>(m_AlnLinksParams[idString].subjRange->GetFrom() + 1,m_AlnLinksParams[idString].subjRange->GetTo() + 1) :
								CRange<TSeqPos>(0,0);					
                        bool flip = (m_AlnLinksParams.count(idString) > 0) ? m_AlnLinksParams[idString].flip : false;
                        int linkout = m_UseLinkoutDB
                            ? CLinkoutDB::GetInstance().GetLinkout(m_AV->GetSeqId(row))
                            : CAlignFormatUtil::GetLinkout(m_AV->GetBioseqHandle(row),m_AV->GetSeqId(row));
                        CAlignFormatUtil::SSeqURLInfo seqUrlInfo("",m_BlastType,m_IsDbNa,m_DbName,m_Rid,
                                                 m_QueryNumber,
                                                 gi,
						                         alnRoInfo->seqidArray[row], 
                                                 linkout,
                                                 m_cur_align,
                                                 true,
                                                 (m_AlignOption & eNewTargetWindow) ? true : false,
                                                 range,
                                                 flip,
                                                 alnRoInfo->taxid[row],
                                                 (m_AlignOption & eShowInfoOnMouseOverSeqid) ? true : false);
                        
                        urlLink = x_GetUrl(&seqUrlInfo, m_AV->GetBioseqHandle(row).GetBioseqCore()->GetId());                        

                        
                        out << urlLink;            
                    }        
                }
                
                bool has_mismatch = false;
                //change the alignment line to identity style
                if (row>0 && m_AlignOption & eShowIdentity){
                    for (int index = j; index < j + (int)actualLineLen && 
                             index < (int)alnRoInfo->sequence[row].size(); index ++){
                        if (alnRoInfo->sequence[row][index] == alnRoInfo->sequence[0][index] &&
                            isalpha((unsigned char) alnRoInfo->sequence[row][index])) {
                            alnRoInfo->sequence[row][index] = k_IdentityChar;           
                        } else if (!has_mismatch) {
                            has_mismatch = true;
                        }        
                    }
                }
                
                //highlight the seqid for pairwise-with-identity format
                if(row>0 && m_AlignOption&eHtml && !(m_AlignOption&eMergeAlign)
                   && m_AlignOption&eShowIdentity && has_mismatch && 
                   (m_AlignOption & eColorDifferentBases)){
                    out<< "<font color = \""<<k_ColorRed<<"\"><b>";         
                }
                out<<alnRoInfo->seqidArray[row]; 
               
                //highlight the seqid for pairwise-with-identity format
                if(row>0 && m_AlignOption&eHtml && !(m_AlignOption&eMergeAlign)
                   && m_AlignOption&eShowIdentity && has_mismatch){
                    out<< "</b></font>";         
                } 
               
                if(urlLink != NcbiEmptyString){
                    //mouse over seqid defline
                    if(m_AlignOption&eHtml &&
                       m_AlignOption&eShowInfoOnMouseOverSeqid) {
                        out << "<span>" <<
                            GetTitle(m_AV->GetBioseqHandle(row)) << "</span>";
                    }
                    out<<"</a>";   
                }
                
                //print out sequence line
                //adjust space between id and start
                CAlignFormatUtil::AddSpace(out, 
                                           alnRoInfo->maxIdLen-alnRoInfo->seqidArray[row].size()+
                                           k_IdStartMargin);
                //not to display start and stop number for empty row
                if (j > 0 && end == prev_stop[row] 
                    || j == 0 && start == 1 && end == 1) {
                    startLen = 0;
                } else {
                    out << start;
                    startLen=NStr::IntToString(start).size();
                }

                CAlignFormatUtil::AddSpace(out, alnRoInfo->maxStartLen-startLen+
                                           k_StartSequenceMargin);
                x_OutputSeq(alnRoInfo->sequence[row], m_AV->GetSeqId(row), j, 
                            (int)actualLineLen, alnRoInfo->frame[row], row,
                            (row > 0 && alnRoInfo->colorMismatch)?true:false,  
                            alnRoInfo->masked_regions[row], out);
                CAlignFormatUtil::AddSpace(out, k_SeqStopMargin);

                 //not to display stop number for empty row in the middle
                if (!(j > 0 && end == prev_stop[row])
                    && !(j == 0 && start == 1 && end == 1)) {
                    out << end;
                }
                
                out<<"\n";
                if(m_AlignOption & eMasterAnchored){//inserts for anchored view
                    bool insertAlready = false;
                    for(list<string>::iterator iter = inserts.begin(); 
                        iter != inserts.end(); iter ++){   
                        if(!insertAlready){
                            if((m_AlignOption&eHtml)
                               &&(m_AlignOption&eMergeAlign) 
                               && (m_AlignOption&eSequenceRetrieval 
                                   && m_CanRetrieveSeq)){
                                char checkboxBuf[200];
                                sprintf(checkboxBuf, 
                                        k_UncheckabeCheckbox.c_str(),
                                        m_QueryNumber);
                                out << checkboxBuf;
                            }
                            CAlignFormatUtil::AddSpace(out, 
                                                       alnRoInfo->maxIdLen
                                                       +k_IdStartMargin
                                                       +alnRoInfo->maxStartLen
                                                       +k_StartSequenceMargin);
                            out << insertPosString<<"\n";
                        }
                        if((m_AlignOption&eHtml)
                           &&(m_AlignOption&eMergeAlign) 
                           && (m_AlignOption&eSequenceRetrieval && m_CanRetrieveSeq)){
                            char checkboxBuf[200];
                            sprintf(checkboxBuf, k_UncheckabeCheckbox.c_str(),
                                    m_QueryNumber);
                            out << checkboxBuf;
                        }
                        CAlignFormatUtil::AddSpace(out, alnRoInfo->maxIdLen
                                                   +k_IdStartMargin
                                                   +alnRoInfo->maxStartLen
                                                   +k_StartSequenceMargin);
                        out<<*iter<<"\n";
                        insertAlready = true;
                    }
                } 
                //display subject sequence feature.
                if(row > 0){ 
                    x_PrintFeatures(alnRoInfo->bioseqFeature[row], row, curRange,
                                    j,(int)actualLineLen, alnRoInfo->maxIdLen,
                                    alnRoInfo->maxStartLen, alnRoInfo->max_feature_num, 
                                    master_feat_str, out);
                }
                //display middle line
                if (row == 0 && ((m_AlignOption & eShowMiddleLine)) 
                    && !(m_AlignOption&eMergeAlign)) {
                    CSeq_id no_id;
                    CAlignFormatUtil::
                        AddSpace(out, alnRoInfo->maxIdLen + k_IdStartMargin
                                 + alnRoInfo->maxStartLen + k_StartSequenceMargin);
                    x_OutputSeq(alnRoInfo->middleLine, no_id, j, (int)actualLineLen, 0,
                                row, false, alnRoInfo->masked_regions[row], out);
                    out<<"\n";
                }
                prev_stop[row] = end; 
            }
            if(!alnRoInfo->seqStarts[row].empty()){ //shouldn't need this check
                alnRoInfo->seqStarts[row].pop_front();
            }
            if(!alnRoInfo->seqStops[row].empty()){
                alnRoInfo->seqStops[row].pop_front();
            }
        }
        out<<"\n";
    }//end of displaying rows    
    string formattedString = CNcbiOstrstreamToString(out);    
    return formattedString;
}

void CDisplaySeqalign::x_PrepareIdentityInfo(SAlnInfo* aln_vec_info)
{
    size_t aln_stop=m_AV->GetAlnStop();
    
    aln_vec_info->match = 0;
    aln_vec_info->positive = 0;
    aln_vec_info->gap = 0;
    aln_vec_info->identity = 0;
    x_FillIdentityInfo(aln_vec_info->alnRowInfo->sequence[0], 
                       aln_vec_info->alnRowInfo->sequence[1], 
                       aln_vec_info->match, 
                       aln_vec_info->positive, 
                       aln_vec_info->alnRowInfo->middleLine);
    if(m_AlignOption & eShowBlastInfo){
        aln_vec_info->identity = CAlignFormatUtil::GetPercentMatch(aln_vec_info->match, (int)aln_stop+1);
        if(aln_vec_info->identity >= k_ColorMismatchIdentity && aln_vec_info->identity <100 &&
               (m_AlignOption & eColorDifferentBases)){
            aln_vec_info->alnRowInfo->colorMismatch = true;
        }
        aln_vec_info->gap = x_GetNumGaps();        
    }
}

void CDisplaySeqalign::x_DisplayAlnvec(CNcbiOstream& out)
{ 
    SAlnRowInfo *alnRoInfo = x_PrepareRowData();    

    string alignRows = x_DisplayRowData(alnRoInfo);
    out << alignRows;
    delete alnRoInfo;
}

CRef<CAlnVec> CDisplaySeqalign::x_GetAlnVecForSeqalign(const CSeq_align& align)
{
    
    //make alnvector
    CRef<CAlnVec> avRef;
    CConstRef<CSeq_align> finalAln;
    if (align.GetSegs().Which() == CSeq_align::C_Segs::e_Std) {
        CRef<CSeq_align> densegAln = align.CreateDensegFromStdseg();
        if (m_AlignOption & eTranslateNucToNucAlignment) { 
            finalAln = densegAln->CreateTranslatedDensegFromNADenseg();
        } else {
            finalAln = densegAln;
        }            
    } else if(align.GetSegs().Which() == 
              CSeq_align::C_Segs::e_Denseg){
        if (m_AlignOption & eTranslateNucToNucAlignment) { 
            finalAln = align.CreateTranslatedDensegFromNADenseg();
        } else {
            finalAln = &align;
        }
    } else if(align.GetSegs().Which() == 
              CSeq_align::C_Segs::e_Dendiag){
        CRef<CSeq_align> densegAln = 
            CAlignFormatUtil::CreateDensegFromDendiag(align);
        if (m_AlignOption & eTranslateNucToNucAlignment) { 
            finalAln = densegAln->CreateTranslatedDensegFromNADenseg();
        } else {
            finalAln = densegAln;
        }
    } else {
        NCBI_THROW(CException, eUnknown, 
                   "Seq-align should be Denseg, Stdseg or Dendiag!");
    }
    CRef<CDense_seg> finalDenseg(new CDense_seg);
    const CTypeConstIterator<CDense_seg> ds = ConstBegin(*finalAln);
    if((ds->IsSetStrands() 
        && ds->GetStrands().front()==eNa_strand_minus) 
       && !(ds->IsSetWidths() && ds->GetWidths()[0] == 3)){
        //show plus strand if master is minus for non-translated case
        finalDenseg->Assign(*ds);
        finalDenseg->Reverse();
        avRef = new CAlnVec(*finalDenseg, m_Scope);   
    } else {
        avRef = new CAlnVec(*ds, m_Scope);
    }    
    
    return avRef;
}


void CDisplaySeqalign::x_FeatSetup(CSeq_align_set &actual_aln_list)
{
    //scope for feature fetching
    if(!(m_AlignOption & eMasterAnchored) 
       && (m_AlignOption & eShowCdsFeature || m_AlignOption 
           & eShowGeneFeature)){
        m_FeatObj = CObjectManager::GetInstance();
        CGBDataLoader::RegisterInObjectManager(*m_FeatObj);
        m_featScope = new CScope(*m_FeatObj);  //for seq feature fetch
        string name = CGBDataLoader::GetLoaderNameFromArgs();
        m_featScope->AddDataLoader(name);
    }   
    m_CanRetrieveSeq = x_GetDbType(actual_aln_list) == eDbTypeNotSet ? false : true;
    if(m_AlignOption & eHtml || m_AlignOption & eDynamicFeature){
        //set config file
        m_ConfigFile = new CNcbiIfstream(".ncbirc");
        m_Reg = new CNcbiRegistry(*m_ConfigFile);
        string feat_file = m_Reg->Get("FEATURE_INFO", "FEATURE_FILE");
        string feat_file_index = m_Reg->Get("FEATURE_INFO",
                                            "FEATURE_FILE_INDEX");
        if(feat_file != NcbiEmptyString && feat_file_index != NcbiEmptyString){
            m_DynamicFeature = new CGetFeature(feat_file, feat_file_index);
        }
    }
}

void CDisplaySeqalign::DisplaySeqalign(CNcbiOstream& out)
{   
    CSeq_align_set actual_aln_list;
    CAlignFormatUtil::ExtractSeqalignSetFromDiscSegs(actual_aln_list, 
                                                     *m_SeqalignSetRef);
    if (actual_aln_list.Get().empty()){
        return;
    }
    
    //sets m_featScope, m_CanRetrieveSeq,m_DynamicFeature
    x_FeatSetup(actual_aln_list);    
    
    if(m_AlignOption & eHtml){  
        out<<"<script src=\"blastResult.js\"></script>";
    }
    //get sequence
    if(m_AlignOption&eSequenceRetrieval && m_AlignOption&eHtml && m_CanRetrieveSeq){ 
        if(m_Ctx->GetRequestValue("NEW_DESIGN").GetValue()!="on")
            out<<s_GetSeqForm((char*)"submitterTop", m_IsDbNa, m_QueryNumber, 
                          x_GetDbType(actual_aln_list),m_DbName, m_Rid.c_str(),
                          s_GetQueryIDFromSeqAlign(actual_aln_list).c_str(),
                          ((m_AlignOption & eDisplayTreeView) ? true: false));
        out<<"<form name=\"getSeqAlignment"<<m_QueryNumber<<"\">\n";
    }
    //begin to display
    int num_align = 0;
    m_cur_align = 0;
    m_currAlignHsp = 0;    
    auto_ptr<CObjectOStream> out2(CObjectOStream::Open(eSerial_AsnText, out));
    //*out2 << *m_SeqalignSetRef;
    //get segs first and get hsp number - m_segs,m_Hsp,m_subjRange
    x_PreProcessSeqAlign(actual_aln_list);
    if(!(m_AlignOption&eMergeAlign)){
        /*pairwise alignment. Note we can't just show each alnment as we go
          because we will need seg information form all hsp's with the same id
          for genome url link.  As a result we show hsp's with the same id 
          as a group*/
     
        CConstRef<CSeq_id> previousId, subid;
        for (CSeq_align_set::Tdata::const_iterator 
                 iter =  actual_aln_list.Get().begin(); 
             iter != actual_aln_list.Get().end() 
                 && num_align<m_NumAlignToShow; iter++, num_align++) {

            //make alnvector
            CRef<CAlnVec> avRef = x_GetAlnVecForSeqalign(**iter);
            
            if(!(avRef.Empty())){
                //Note: do not switch the set order per calnvec specs.
                avRef->SetGenCode(m_SlaveGeneticCode);
                avRef->SetGenCode(m_MasterGeneticCode, 0);
                try{
                    const CBioseq_Handle& handle = avRef->GetBioseqHandle(1);
                    if(handle){
                      
                        //save the current alnment regardless
                        CRef<SAlnInfo> alnvecInfo(new SAlnInfo);
                        int num_ident;
                        CAlignFormatUtil::GetAlnScores(**iter, 
                                                       alnvecInfo->score, 
                                                       alnvecInfo->bits, 
                                                       alnvecInfo->evalue, 
                                                       alnvecInfo->sum_n, 
                                                       num_ident,
                                                       alnvecInfo->use_this_gi,
                                                       alnvecInfo->comp_adj_method);
                        alnvecInfo->alnvec = avRef;
                       
                        subid=&(avRef->GetSeqId(1));
                        if(!previousId.Empty() && 
                           !subid->Match(*previousId)){
                            m_Scope.RemoveFromHistory(m_Scope.
                                                      GetBioseqHandle(*previousId));
                                                      //release memory 
                        }
                        
                        x_DisplayAlnvecInfo(out, alnvecInfo,
                                            previousId.Empty() || 
                                            !subid->Match(*previousId));
                       
                        previousId = subid;
                    }                
                } catch (const CException&){
                    out << "Sequence with id "
                        << (avRef->GetSeqId(1)).GetSeqIdString().c_str() 
                        <<" no longer exists in database...alignment skipped\n";
                    continue;
                }
            }
        } 

    } else if(m_AlignOption&eMergeAlign){ //multiple alignment
        vector< CRef<CAlnMix> > mix(k_NumFrame); 
        //each for one frame for translated alignment
        for(int i = 0; i < k_NumFrame; i++){
            mix[i] = new CAlnMix(m_Scope);
        }        
        num_align = 0;
        vector<CRef<CSeq_align_set> > alnVector(k_NumFrame);
        for(int i = 0; i <  k_NumFrame; i ++){
            alnVector[i] = new CSeq_align_set;
        }
        for (CSeq_align_set::Tdata::const_iterator 
                 alnIter = actual_aln_list.Get().begin(); 
             alnIter != actual_aln_list.Get().end() 
                 && num_align<m_NumAlignToShow; alnIter ++, num_align++) {

            const CBioseq_Handle& subj_handle = 
                m_Scope.GetBioseqHandle((*alnIter)->GetSeq_id(1));
            if(subj_handle){
                //need to convert to denseg for stdseg
                if((*alnIter)->GetSegs().Which() == CSeq_align::C_Segs::e_Std) {
                    CTypeConstIterator<CStd_seg> ss = ConstBegin(**alnIter); 
                    CRef<CSeq_align> convertedDs = 
                        (*alnIter)->CreateDensegFromStdseg();
                    if((convertedDs->GetSegs().GetDenseg().IsSetWidths() 
                        && convertedDs->GetSegs().GetDenseg().GetWidths()[0] == 3)
                       || m_AlignOption & eTranslateNucToNucAlignment){
                        //only do this for translated master
                        int frame = s_GetStdsegMasterFrame(*ss, m_Scope);
                        switch(frame){
                        case 1:
                            alnVector[0]->Set().push_back(convertedDs);
                            break;
                        case 2:
                            alnVector[1]->Set().push_back(convertedDs);
                            break;
                        case 3:
                            alnVector[2]->Set().push_back(convertedDs);
                            break;
                        case -1:
                            alnVector[3]->Set().push_back(convertedDs);
                            break;
                        case -2:
                            alnVector[4]->Set().push_back(convertedDs);
                            break;
                        case -3:
                            alnVector[5]->Set().push_back(convertedDs);
                            break;
                        default:
                            break;
                        }
                    }
                    else {
                        alnVector[0]->Set().push_back(convertedDs);
                    }
                } else if((*alnIter)->GetSegs().Which() == CSeq_align::C_Segs::
                          e_Denseg){
                    alnVector[0]->Set().push_back(*alnIter);
                } else if((*alnIter)->GetSegs().Which() == CSeq_align::C_Segs::
                          e_Dendiag){
                    alnVector[0]->Set().\
                        push_back(CAlignFormatUtil::CreateDensegFromDendiag(**alnIter));
                } else {
                    NCBI_THROW(CException, eUnknown, 
                               "Input Seq-align should be Denseg, Stdseg or Dendiag!");
                }
            }
        }
        for(int i = 0; i < (int)alnVector.size(); i ++){
            bool hasAln = false;
            for(CTypeConstIterator<CSeq_align> 
                    alnRef = ConstBegin(*alnVector[i]); alnRef; ++alnRef){
                CTypeConstIterator<CDense_seg> ds = ConstBegin(*alnRef);
                //*out2 << *ds;      
                try{
                    if (m_AlignOption & eTranslateNucToNucAlignment) {         
                        mix[i]->Add(*ds, CAlnMix::fForceTranslation);
                    } else {
                        if (ds->IsSetWidths() &&
                            ds->GetWidths()[0] == 3 && 
                            ds->IsSetStrands() && 
                            ds->GetStrands().front()==eNa_strand_minus){
                            mix[i]->Add(*ds, CAlnMix::fNegativeStrand);
                        } else {
                            mix[i]->Add(*ds);
                        }
                    }
                } catch (const CException& e){
                    out << "Warning: " << e.what();
                    continue;
                }
                hasAln = true;
            }
            if(hasAln){
                //    *out2<<*alnVector[i];
                mix[i]->Merge(CAlnMix::fMinGap 
                              | CAlnMix::fQuerySeqMergeOnly 
                              | CAlnMix::fFillUnalignedRegions);  
                //*out2<<mix[i]->GetDenseg();
            }
        }
        
        int numDistinctFrames = 0;
        for(int i = 0; i < (int)alnVector.size(); i ++){
            if(!alnVector[i]->Get().empty()){
                numDistinctFrames ++;
            }
        }
        out<<"\n";
        for(int i = 0; i < k_NumFrame; i ++){
            try{
                CRef<CAlnVec> avRef (new CAlnVec (mix[i]->GetDenseg(), 
                                                  m_Scope));
                avRef->SetGenCode(m_SlaveGeneticCode);
                avRef->SetGenCode(m_MasterGeneticCode, 0);
                m_AV = avRef;
                
                if(numDistinctFrames > 1){
                    out << "For reading frame " << k_FrameConversion[i] 
                        << " of query sequence:\n\n";
                }
                x_DisplayAlnvec(out);
            } catch (CException e){
                continue;
            }
        } 
    }
    if(m_AlignOption&eSequenceRetrieval && m_AlignOption&eHtml && m_CanRetrieveSeq){
        out<<"</form>\n";        
        if(m_Ctx->GetRequestValue("NEW_DESIGN").GetValue()!="on")
            out<<s_GetSeqForm((char*)"submitterBottom", m_IsDbNa,
                              m_QueryNumber, x_GetDbType(actual_aln_list),
                              m_DbName, m_Rid.c_str(),
                              s_GetQueryIDFromSeqAlign(actual_aln_list).c_str(),
                          ((m_AlignOption & eDisplayTreeView) ? true: false));
    }
}


void CDisplaySeqalign::x_FillIdentityInfo(const string& sequence_standard,
                                          const string& sequence , 
                                          int& match, int& positive, 
                                          string& middle_line) 
{
    match = 0;
    positive = 0;
    int min_length=min<int>((int)sequence_standard.size(), (int)sequence.size());
    if(m_AlignOption & eShowMiddleLine){
        middle_line = sequence;
    }
    for(int i=0; i<min_length; i++){
        if(sequence_standard[i]==sequence[i]){
            if(m_AlignOption & eShowMiddleLine){
                if(m_MidLineStyle == eBar ) {
                    middle_line[i] = '|';
                } else if (m_MidLineStyle == eChar){
                    middle_line[i] = sequence[i];
                }
            }
            match ++;
        } else {
            if ((m_AlignType&eProt) 
                && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){  
                positive ++;
                if(m_AlignOption & eShowMiddleLine){
                    if (m_MidLineStyle == eChar){
                        middle_line[i] = '+';
                    }
                }
            } else {
                if (m_AlignOption & eShowMiddleLine){
                    middle_line[i] = ' ';
                }
            }    
        }
    }  
}

CDisplaySeqalign::SAlnDispParams *CDisplaySeqalign::x_FillAlnDispParams(const CRef< CBlast_def_line > &bdl,
								   list<int>& use_this_gi,
								   int firstGi,
								   bool isNa,
								  int seqLength)
{
    SAlnDispParams *alnDispParams = NULL;

	const list<CRef<CSeq_id> > ids = bdl->GetSeqid();
	int gi =  s_GetGiForSeqIdList(ids);
    int gi_in_use_this_gi = 0;
    
    ITERATE(list<int>, iter_gi, use_this_gi){
        if(gi == *iter_gi){
            gi_in_use_this_gi = *iter_gi;
            break;
        }
    }
	if(use_this_gi.empty() || gi_in_use_this_gi > 0) {
		alnDispParams = new SAlnDispParams();
		alnDispParams->gi =  gi;    		
		alnDispParams->seqID = FindBestChoice(ids, CSeq_id::WorstRank);		
		alnDispParams->label =  CAlignFormatUtil::GetLabel(alnDispParams->seqID);//Just accession without db part like ref| or pdbd|
		if(m_AlignOption&eHtml){
			int taxid = 0;
			string type_temp = m_BlastType;
			type_temp = NStr::TruncateSpaces(NStr::ToLower(type_temp));
			if((type_temp == "mapview" || 
                            type_temp == "mapview_prev" || 
                            type_temp == "gsfasta" || type_temp == "gsfasta_prev") && 
                            bdl->IsSetTaxid() && 
                            bdl->CanGetTaxid()){
				taxid = bdl->GetTaxid();
			}            		
			string idString = m_AV->GetSeqId(1).GetSeqIdString();			
			CRange<TSeqPos> range = (m_AlnLinksParams.count(idString) > 0 && m_AlnLinksParams[idString].subjRange) ?
						CRange<TSeqPos>(m_AlnLinksParams[idString].subjRange->GetFrom() + 1,m_AlnLinksParams[idString].subjRange->GetTo() + 1) :
						CRange<TSeqPos>(0,0);					
			bool flip = (m_AlnLinksParams.count(idString) > 0) ? m_AlnLinksParams[idString].flip : false;
            int linkout = m_UseLinkoutDB
                ? CLinkoutDB::GetInstance().GetLinkout(gi) 
                : CAlignFormatUtil::GetLinkout(*bdl);
			CAlignFormatUtil::SSeqURLInfo seqUrlInfo("",m_BlastType,m_IsDbNa,m_DbName,m_Rid,
                                                 m_QueryNumber,
                                                 gi_in_use_this_gi, 
						                         alnDispParams->label,
                                                 linkout,
                                                 m_cur_align,
                                                 true,
                                                 (m_AlignOption & eNewTargetWindow) ? true : false,
                                                 range,
                                                 flip,
                                                 taxid,
                                                 (m_AlignOption & eShowInfoOnMouseOverSeqid) ? true : false);
                        
			alnDispParams->id_url =  x_GetUrl(&seqUrlInfo, ids,m_AlignTemplates != NULL);                                                
		}
		
		if(m_AlignOption&eLinkout){
		
			int linkout = m_UseLinkoutDB
                ? CLinkoutDB::GetInstance().GetLinkout(gi)
                : CAlignFormatUtil::GetLinkout((*bdl));
			string user_url = m_Reg->Get(m_BlastType,"TOOL_URL");
			list<string> linkout_url =  CAlignFormatUtil::
                                GetLinkoutUrl(linkout, ids,
                                              m_Rid,
                                              m_CddRid, m_EntrezTerm,
                                              isNa, 
                                              firstGi,
                                              false, true, m_cur_align);                            
			ITERATE(list<string>, iter_linkout, linkout_url){
				alnDispParams->linkoutStr += *iter_linkout;
			}
			if(seqLength > k_GetSubseqThreshhold){
				alnDispParams->dumpGnlUrl = x_GetDumpgnlLink(ids);                                
			}
        
		}
		if(bdl->IsSetTitle()){
			alnDispParams->title = bdl->GetTitle();
		}
	}
	return alnDispParams;
}



CDisplaySeqalign::SAlnDispParams *CDisplaySeqalign::x_FillAlnDispParams(const CBioseq_Handle& bsp_handle)
{
    SAlnDispParams *alnDispParams = new SAlnDispParams();
	alnDispParams->gi = FindGi(bsp_handle.GetBioseqCore()->GetId());
	alnDispParams->seqID = FindBestChoice(bsp_handle.GetBioseqCore()->GetId(),CSeq_id::WorstRank);
	alnDispParams->label =  CAlignFormatUtil::GetLabel(alnDispParams->seqID);
	if(m_AlignOption&eHtml){            		
			string idString = m_AV->GetSeqId(1).GetSeqIdString();			
			CRange<TSeqPos> range = (m_AlnLinksParams.count(idString) > 0 && m_AlnLinksParams[idString].subjRange) ?
						CRange<TSeqPos>(m_AlnLinksParams[idString].subjRange->GetFrom() + 1,m_AlnLinksParams[idString].subjRange->GetTo() + 1) :
						CRange<TSeqPos>(0,0);					
			bool flip = (m_AlnLinksParams.count(idString) > 0) ? m_AlnLinksParams[idString].flip : false;
			CAlignFormatUtil::SSeqURLInfo seqUrlInfo("",m_BlastType,m_IsDbNa,m_DbName,m_Rid,
                                                 m_QueryNumber,
                                                 alnDispParams->gi, 
						                         alnDispParams->label,
                                                 0,
                                                 m_cur_align,
                                                 true,
                                                 (m_AlignOption & eNewTargetWindow) ? true : false,
                                                 range,
                                                 flip,
                                                 0,
                                                 (m_AlignOption & eShowInfoOnMouseOverSeqid) ? true : false);
                        
			alnDispParams->id_url = x_GetUrl(&seqUrlInfo, bsp_handle.GetBioseqCore()->GetId(),m_AlignTemplates != NULL);
	}
			
	alnDispParams->title = GetTitle(bsp_handle);		
	
	return alnDispParams;
}

string
CDisplaySeqalign::x_PrintDefLine(const CBioseq_Handle& bsp_handle,SAlnInfo* aln_vec_info)                                 
{
    CNcbiOstrstream out;                
    /* Facilitates comparing formatted output using diff */
    static string kLengthString("Length=");
#ifdef CTOOLKIT_COMPATIBLE
    static bool value_set = false;
    if ( !value_set ) {
        if (getenv("CTOOLKIT_COMPATIBLE")) {
            kLengthString.assign("          Length = ");
        }
        value_set = true;
    }
#endif /* CTOOLKIT_COMPATIBLE */
  
    if(bsp_handle){
        const CRef<CSeq_id> wid =
            FindBestChoice(bsp_handle.GetBioseqCore()->GetId(), 
                           CSeq_id::WorstRank);
    
        const CRef<CBlast_def_line_set> bdlRef 
            =  CAlignFormatUtil::GetBlastDefline(bsp_handle);
        const list< CRef< CBlast_def_line > >& bdl = bdlRef->Get();
        bool isFirst = true;
        int firstGi = 0;

        m_cur_align++;
    
        if(bdl.empty()){ //no blast defline struct, should be no such case now
            //actually not so fast...as we now fetch from entrez even when it's not in blast db
            //there is no blast defline in such case.
			SAlnDispParams *alnDispParams = x_FillAlnDispParams(bsp_handle);
            out << ">";
            if ((m_AlignOption&eSequenceRetrieval)
                && (m_AlignOption&eHtml) && m_CanRetrieveSeq && isFirst) {
                char buf[512];
                sprintf(buf, k_Checkbox.c_str(), alnDispParams->gi > 0 ?
                        NStr::IntToString(alnDispParams->gi).c_str() : alnDispParams->label.c_str(),
                        m_QueryNumber);
                out << buf;
            }
                
            if(m_AlignOption&eHtml){               
		        
                aln_vec_info->id_label = (alnDispParams->gi != 0) ? NStr::IntToString(alnDispParams->gi) : alnDispParams->label;                      

                out<<alnDispParams->id_url;
            }
                
            if(m_AlignOption&eShowGi && alnDispParams->gi > 0){
                out<<"gi|"<<alnDispParams->gi<<"|";
                    }     
            if(!(alnDispParams->seqID->AsFastaString().find("gnl|BL_ORD_ID") 
                 != string::npos)){
                alnDispParams->seqID->WriteAsFasta(out);
            }
            if(m_AlignOption&eHtml){
                if(alnDispParams->id_url != NcbiEmptyString){
                    out<<"</a>";
                }
                if(alnDispParams->gi != 0){
                    out<<"<a name="<<alnDispParams->gi<<"></a>";                    
                } else {
                    out<<"<a name="<<alnDispParams->seqID->GetSeqIdString()<<"></a>";                    
                }
            }
            out <<" ";
            s_WrapOutputLine(out, (m_AlignOption&eHtml) ? 
                             CHTMLHelper::HTMLEncode(alnDispParams->title) :
                             alnDispParams->title);     
                
            out<<"\n";
            delete alnDispParams;    
        } else {
            //print each defline 
            bool bMultipleDeflines = false;
            int numBdl = 0;
            for(list< CRef< CBlast_def_line > >::const_iterator 
                    iter = bdl.begin(); iter != bdl.end(); iter++){
                numBdl++;
				SAlnDispParams *alnDispParams = x_FillAlnDispParams(*iter,
																	aln_vec_info->use_this_gi,
																	firstGi,
																	bsp_handle.GetBioseqCore()->IsNa(), 
																	(int)bsp_handle.GetBioseqLength());                                                                    

                if(alnDispParams) {
                
                    if(isFirst){
                        out << ">";                  
                    } else{
                        out << " ";
                        if (m_AlignOption&eHtml && (int)(bdl.size()) > k_MaxDeflinesToShow && numBdl == k_MinDeflinesToShow + 1){ 
                            //Show first 3 deflines out of 8 or more, hide the rest
                            string mdlTag = aln_vec_info->id_label;
                            //string mdlTag = id_label  + "_" + NStr::IntToString(m_cur_align);                        
                            out << "<a href=\"#\" title=\"Other sequence titles\"  onmouseover=\"showInfo(this)\" class=\"resArrowLinkW mdl hiding\" id=\"" <<
                                mdlTag << "\">" << bdl.size() - k_MinDeflinesToShow << " more sequence titles" << "</a>\n";
                        
                            out << " <div id=\"" << "info_" << mdlTag << "\" class=\"helpbox mdlbox hidden\">";
                            bMultipleDeflines = true;
                        }
                    }
                    
                    if(isFirst){
                        firstGi = alnDispParams->gi;
                    }
                    if ((m_AlignOption&eSequenceRetrieval)
                        && (m_AlignOption&eHtml) && m_CanRetrieveSeq && isFirst) {
                        char buf[512];
                        sprintf(buf, k_Checkbox.c_str(), alnDispParams->gi > 0 ?
                                NStr::IntToString(alnDispParams->gi).c_str() : alnDispParams->label.c_str(),
                                m_QueryNumber);
                                out << buf;
                    }
                
                    if(m_AlignOption&eHtml){
                        out<< alnDispParams->id_url;
                    }
                
                    if(m_AlignOption&eShowGi && alnDispParams->gi > 0){
                        out<<"gi|"<<alnDispParams->gi<<"|";
                    }     
                    if(!(alnDispParams->seqID->AsFastaString().find("gnl|BL_ORD_ID") 
                         != string::npos)){
                        alnDispParams->seqID->WriteAsFasta(out);
                    }
                    if(m_AlignOption&eHtml){
                        if(alnDispParams->id_url != NcbiEmptyString){
                            out<<"</a>";
                        }
                        if(alnDispParams->gi != 0){
                            out<<"<a name="<<alnDispParams->gi<<"></a>";
                            aln_vec_info->id_label = NStr::IntToString(alnDispParams->gi);
                        } else {
                            out<<"<a name="<<alnDispParams->seqID->GetSeqIdString()<<"></a>";
                            aln_vec_info->id_label = alnDispParams->label;
                        }
                        if(m_AlignOption&eLinkout){
                            
                            out <<" ";
                            out << alnDispParams->linkoutStr;
							if(!alnDispParams->dumpGnlUrl.empty()) {                           
                            
                                out<<alnDispParams->dumpGnlUrl;
                            }
						}
                    }
                
                    out <<" ";
					if(!alnDispParams->title.empty()) {                    
                        s_WrapOutputLine(out, (m_AlignOption&eHtml) ? 
                                         CHTMLHelper::
                                         HTMLEncode(alnDispParams->title) :
                                         alnDispParams->title);     
                    }
                    out<<"\n";
                    isFirst = false;
                    delete alnDispParams;
                }
            }
            if(m_AlignOption&eHtml && bMultipleDeflines) {
                out << "</div>";
            }                      
        }
    }      
    out<<kLengthString<<bsp_handle.GetBioseqLength()<<"\n";    
    string formattedString = CNcbiOstrstreamToString(out);    
    return formattedString;
}


void CDisplaySeqalign::x_OutputSeq(string& sequence, const CSeq_id& id, 
                                   int start, int len, int frame, int row,
                                   bool color_mismatch,
                                   const TSAlnSeqlocInfoList& loc_list, 
                                   CNcbiOstream& out) const 
{
    _ASSERT((int)sequence.size() > start);
    list<CRange<int> > actualSeqloc;
    string actualSeq = sequence.substr(start, len);
    
    if(id.Which() != CSeq_id::e_not_set){
        /*only do this for sequence but not for others like middle line,
          features*/
        ITERATE(TSAlnSeqlocInfoList, iter, loc_list) {
            int from=(*iter)->aln_range.GetFrom();
            int to=(*iter)->aln_range.GetTo();
            int locFrame = (*iter)->seqloc->GetFrame();
            if(id.Match((*iter)->seqloc->GetInterval().GetId()) 
               && locFrame == frame){
                bool isFirstChar = true;
                CRange<int> eachSeqloc(0, 0);
                //go through each residule and mask it
                for (int i=max<int>(from, start); 
                     i<=min<int>(to, start+len); i++){
                    //store seqloc start for font tag below
                    if ((m_AlignOption & eHtml) && isFirstChar){         
                        isFirstChar = false;
                        eachSeqloc.Set(i, eachSeqloc.GetTo());
                    }
                    if (m_SeqLocChar==eX){
                        if(isalpha((unsigned char) actualSeq[i-start])){
                            actualSeq[i-start]='X';
                        }
                    } else if (m_SeqLocChar==eN){
                        actualSeq[i-start]='n';
                    } else if (m_SeqLocChar==eLowerCase){
                        actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
                    }
                    //store seqloc start for font tag below
                    if ((m_AlignOption & eHtml) 
                        && i == min<int>(to, start+len)){ 
                        eachSeqloc.Set(eachSeqloc.GetFrom(), i);
                    }
                }
                if(!(eachSeqloc.GetFrom()==0&&eachSeqloc.GetTo()==0)){
                    actualSeqloc.push_back(eachSeqloc);
                }
            }
        }
    }
    
    if(actualSeqloc.empty()){//no need to add font tag
        if((m_AlignOption & eColorDifferentBases) && (m_AlignOption & eHtml)
           && color_mismatch && (m_AlignOption & eShowIdentity)){
            //color the mismatches. Only for rows without mask. 
            //Otherwise it may confilicts with mask font tag.
            s_ColorDifferentBases(actualSeq, k_IdentityChar, out);
        } else {
            out<<actualSeq;
        }
    } else {//now deal with font tag for mask for html display    
        bool endTag = false;
        bool numFrontTag = 0;
        for (int i = 0; i < (int)actualSeq.size(); i ++){
            for (list<CRange<int> >::iterator iter=actualSeqloc.begin(); 
                 iter!=actualSeqloc.end(); iter++){
                int from = (*iter).GetFrom() - start;
                int to = (*iter).GetTo() - start;
                //start tag
                if(from == i){
                    out<<"<font color=\""+color[m_SeqLocColor]+"\">";
                    numFrontTag = 1;
                }
                //need to close tag at the end of mask or end of sequence
                if(to == i || i == (int)actualSeq.size() - 1 ){
                    endTag = true;
                }
            }
            out<<actualSeq[i];
            if(endTag && numFrontTag == 1){
                out<<"</font>";
                endTag = false;
                numFrontTag = 0;
            }
        }
    }
}


int CDisplaySeqalign::x_GetNumGaps() 
{
    int gap = 0;
    for (int row=0; row<m_AV->GetNumRows(); row++) {
        CRef<CAlnMap::CAlnChunkVec> chunk_vec 
            = m_AV->GetAlnChunks(row, m_AV->GetSeqAlnRange(0));
        for (int i=0; i<chunk_vec->size(); i++) {
            CConstRef<CAlnMap::CAlnChunk> chunk = (*chunk_vec)[i];
            if (chunk->IsGap()) {
                gap += (chunk->GetAlnRange().GetTo() 
                        - chunk->GetAlnRange().GetFrom() + 1);
            }
        }
    }
    return gap;
}


void CDisplaySeqalign::x_GetFeatureInfo(TSAlnFeatureInfoList& feature,
                                        CScope& scope, 
                                        CSeqFeatData::E_Choice choice,
                                        int row, string& sequence,
                                        list<list<CRange<TSeqPos> > >& feat_range_list,
                                        list<ENa_strand>& feat_seq_strand,
                                        bool fill_feat_range ) const 
{
    //Only fetch features for seq that has a gi unless it's master seq
    const CSeq_id& id = m_AV->GetSeqId(row);
    
    int gi_temp = FindGi(m_AV->GetBioseqHandle(row).GetBioseqCore()->GetId());
    if(gi_temp > 0 || row == 0){
        const CBioseq_Handle& handle = scope.GetBioseqHandle(id);
        if(handle){
            TSeqPos seq_start = m_AV->GetSeqPosFromAlnPos(row, 0);
            TSeqPos seq_stop = m_AV->GetSeqPosFromAlnPos(row, m_AV->GetAlnStop());
            CRef<CSeq_loc> loc_ref =
                handle.
                GetRangeSeq_loc(min(seq_start, seq_stop),
                                max(seq_start, seq_stop));
            for (CFeat_CI feat(scope, *loc_ref, choice); feat; ++feat) {
                const CSeq_loc& loc = feat->GetLocation();
                bool has_id = false;
                list<CSeq_loc_CI::TRange> isolated_range;
                ENa_strand feat_strand = eNa_strand_plus, prev_strand = eNa_strand_plus;
                bool first_loc = true, mixed_strand = false, mix_loc = false;
                CRange<TSeqPos> feat_seq_range;
                TSeqPos other_seqloc_length = 0;
                //isolate the seqloc corresponding to feature
                //as this is easier to manipulate and remove seqloc that is
                //not from the bioseq we are dealing with
                for(CSeq_loc_CI loc_it(loc); loc_it; ++loc_it){
                    const CSeq_id& id_it = loc_it.GetSeq_id();
                    if(IsSameBioseq(id_it, id, &scope)){
                        isolated_range.push_back(loc_it.GetRange());
                        if(first_loc){
                            feat_seq_range = loc_it.GetRange();
                        } else {
                            feat_seq_range += loc_it.GetRange();
                        }
                        has_id = true;
                        if(loc_it.IsSetStrand()){
                            feat_strand = loc_it.GetStrand();
                            if(feat_strand != eNa_strand_plus && 
                               feat_strand != eNa_strand_minus){
                                feat_strand = eNa_strand_plus;
                            }
                        } else {
                            feat_strand = eNa_strand_plus;
                        }
                   
                        if(!first_loc && prev_strand != feat_strand){
                            mixed_strand = true;
                        }
                        first_loc = false;
                        prev_strand = feat_strand;
                    } else {
                        //if seqloc has other seqids then need to remove other 
                        //seqid encoded amino acids in the front later
                        if (first_loc) {
                            other_seqloc_length += loc_it.GetRange().GetLength();
                            mix_loc = true;
                        }
                    }
                }
                //give up if mixed strand or no id
                if(!has_id || mixed_strand){
                    continue;
                }
               
                string featLable = NcbiEmptyString;
                string featId;
                char feat_char = ' ';
                string alternativeFeatStr = NcbiEmptyString;            
                TSeqPos feat_aln_from = 0;
                TSeqPos feat_aln_to = 0;
                TSeqPos actual_feat_seq_start = 0, actual_feat_seq_stop = 0;
                feature::GetLabel(feat->GetOriginalFeature(), &featLable, 
                                  feature::eBoth, &scope);
                featId = featLable.substr(0, k_FeatureIdLen); //default
                TSeqPos aln_stop = m_AV->GetAlnStop();  
                CRef<SAlnFeatureInfo> featInfo;
               
                //find the actual feature sequence start and stop 
                if(m_AV->IsPositiveStrand(row)){
                    actual_feat_seq_start = 
                        max(feat_seq_range.GetFrom(), seq_start);
                    actual_feat_seq_stop = 
                        min(feat_seq_range.GetTo(), seq_stop);
                    
                } else {
                    actual_feat_seq_start = 
                        min(feat_seq_range.GetTo(), seq_start);
                    actual_feat_seq_stop =
                        max(feat_seq_range.GetFrom(), seq_stop);
                }
                //the feature alignment positions
                feat_aln_from = 
                    m_AV->GetAlnPosFromSeqPos(row, actual_feat_seq_start);
                feat_aln_to = 
                    m_AV->GetAlnPosFromSeqPos(row, actual_feat_seq_stop);
                if(choice == CSeqFeatData::e_Gene){
                    featInfo.Reset(new SAlnFeatureInfo); 
                    feat_char = '^';
                    
                } else if(choice == CSeqFeatData::e_Cdregion){
                     
                    string raw_cdr_product = 
                        s_GetCdsSequence(m_SlaveGeneticCode, feat, scope,
                                         isolated_range, handle, feat_strand,
                                         featId, other_seqloc_length%3 == 0 ?
                                         0 : 3 - other_seqloc_length%3,
                                         mix_loc);
                    if(raw_cdr_product == NcbiEmptyString){
                        continue;
                    }
                    featInfo.Reset(new SAlnFeatureInfo);
                          
                    //line represents the amino acid line starting covering 
                    //the whole alignment.  The idea is if there is no feature
                    //in some range, then fill it with space and this won't
                    //be shown 
                    
                    string line(aln_stop+1, ' '); 
                    //pre-fill all cds region with intron char
                    for (TSeqPos i = feat_aln_from; i <= feat_aln_to; i ++){
                        line[i] = k_IntronChar;
                    }
                    
                    //get total coding length
                    TSeqPos total_coding_len = 0;
                    ITERATE(list<CSeq_loc_CI::TRange>, iter, isolated_range){
                        total_coding_len += iter->GetLength(); 
                    }
                  
                    //fill concatenated exon (excluding intron)
                    //with product
                    //this is will be later used to
                    //fill the feature line
                    char gap_char = m_AV->GetGapChar(row);
                    string concat_exon = 
                        s_GetConcatenatedExon(feat, feat_strand, 
                                              isolated_range,
                                              total_coding_len,
                                              raw_cdr_product,
                                              other_seqloc_length%3 == 0 ?
                                              0 : 3 - other_seqloc_length%3);
                    
                   
                    //fill slave feature info to make putative feature for
                    //master sequence
                    if (fill_feat_range) {
                        list<CRange<TSeqPos> > master_feat_range;
                        ENa_strand master_strand = eNa_strand_plus;
                        s_MapSlaveFeatureToMaster(master_feat_range, master_strand,
                                                  feat, isolated_range,
                                                  feat_strand, m_AV, row,
                                                  other_seqloc_length%3 == 0 ?
                                                  0 : 
                                                  3 - other_seqloc_length%3);
                        if(!(master_feat_range.empty())) {
                            feat_range_list.push_back(master_feat_range); 
                            feat_seq_strand.push_back(master_strand);
                        } 
                    }
                       
                    
                    TSeqPos feat_aln_start_totalexon = 0;
                    TSeqPos prev_feat_aln_start_totalexon = 0;
                    TSeqPos prev_feat_seq_stop = 0;
                    TSeqPos intron_size = 0;
                    bool is_first = true;
                    bool  is_first_exon_start = true;
               
                    //here things get complicated a bit. The idea is fill the
                    //whole feature line in alignment coordinates with
                    //amino acid on the second base of a condon

                    //go through the feature seqloc and fill the feature line
                    
                    //Need to reverse the seqloc order for minus strand
                    if(feat_strand == eNa_strand_minus){
                        isolated_range.reverse(); 
                    }
                    
                    ITERATE(list<CSeq_loc_CI::TRange>, iter, isolated_range){
                        //intron refers to the distance between two exons
                        //i.e. each seqloc is an exon
                        //intron needs to be skipped
                        if(!is_first){
                            intron_size += iter->GetFrom() 
                                - prev_feat_seq_stop - 1;
                        }
                        CRange<TSeqPos> actual_feat_seq_range =
                            loc_ref->GetTotalRange().
                            IntersectionWith(*iter);          
                        if(!actual_feat_seq_range.Empty()){
                            //the sequence start position in aln coordinates
                            //that has a feature
                            TSeqPos feat_aln_start;
                            TSeqPos feat_aln_stop;
                            if(m_AV->IsPositiveStrand(row)){
                                feat_aln_start = 
                                    m_AV->
                                    GetAlnPosFromSeqPos
                                    (row, actual_feat_seq_range.GetFrom());
                                feat_aln_stop
                                    = m_AV->GetAlnPosFromSeqPos
                                    (row, actual_feat_seq_range.GetTo());
                            } else {
                                feat_aln_start = 
                                    m_AV->
                                    GetAlnPosFromSeqPos
                                    (row, actual_feat_seq_range.GetTo());
                                feat_aln_stop
                                    = m_AV->GetAlnPosFromSeqPos
                                    (row, actual_feat_seq_range.GetFrom());
                            }
                            //put actual amino acid on feature line
                            //in aln coord 
                            for (TSeqPos i = feat_aln_start;
                                 i <= feat_aln_stop;  i ++){  
                                    if(sequence[i] != gap_char){
                                        //the amino acid position in 
                                        //concatanated exon that corresponds
                                        //to the sequence position
                                        //note intron needs to be skipped
                                        //as it does not have cds feature
                                        TSeqPos product_adj_seq_pos
                                            = m_AV->GetSeqPosFromAlnPos(row, i) - 
                                            intron_size - feat_seq_range.GetFrom();
                                        if(product_adj_seq_pos < 
                                           concat_exon.size()){
                                            //fill the cds feature line with
                                            //actual amino acids
                                            line[i] = 
                                                concat_exon[product_adj_seq_pos];
                                            //get the exon start position
                                            //note minus strand needs to be
                                            //counted backward
                                            if(m_AV->IsPositiveStrand(row)){
                                                //don't count gap 
                                                if(is_first_exon_start &&
                                                   isalpha((unsigned char) line[i])){
                                                    if(feat_strand == eNa_strand_minus){ 
                                                        feat_aln_start_totalexon = 
                                                            concat_exon.size()
                                                            - product_adj_seq_pos + 1;
                                                        is_first_exon_start = false;
                                                        
                                                    } else {
                                                        feat_aln_start_totalexon = 
                                                            product_adj_seq_pos;
                                                        is_first_exon_start = false;
                                                    }
                                                }
                                                
                                            } else {
                                                if(feat_strand == eNa_strand_minus){ 
                                                    if(is_first_exon_start && 
                                                       isalpha((unsigned char) line[i])){
                                                        feat_aln_start_totalexon = 
                                                            concat_exon.size()
                                                            - product_adj_seq_pos + 1;
                                                        is_first_exon_start = false;
                                                        prev_feat_aln_start_totalexon =
                                                        feat_aln_start_totalexon;
                                                    }
                                                    if(!is_first_exon_start){
                                                        //need to get the
                                                        //smallest start as
                                                        //seqloc list is
                                                        //reversed
                                                        feat_aln_start_totalexon =
                                                            min(TSeqPos(concat_exon.size()
                                                                        - product_adj_seq_pos + 1), 
                                                                prev_feat_aln_start_totalexon);
                                                        prev_feat_aln_start_totalexon =
                                                        feat_aln_start_totalexon;  
                                                    }
                                                } else {
                                                    feat_aln_start_totalexon = 
                                                        max(prev_feat_aln_start_totalexon,
                                                            product_adj_seq_pos); 
                                                    
                                                    prev_feat_aln_start_totalexon =
                                                        feat_aln_start_totalexon;
                                                }
                                            }
                                        }
                                    } else { //adding gap
                                        line[i] = ' '; 
                                    }                         
                               
                            }                      
                        }
                        
                        prev_feat_seq_stop = iter->GetTo();  
                        is_first = false;
                    }                 
                    alternativeFeatStr = line;
                    s_FillCdsStartPosition(line, concat_exon, m_LineLen,
                                           feat_aln_start_totalexon,
                                           m_AV->IsPositiveStrand(row) ?
                                           eNa_strand_plus : eNa_strand_minus,
                                           feat_strand, featInfo->feature_start); 
                 
                }
                
                if(featInfo){
                    x_SetFeatureInfo(featInfo, *loc_ref,
                                     feat_aln_from, feat_aln_to, aln_stop, 
                                     feat_char, featId, alternativeFeatStr);  
                    feature.push_back(featInfo);
                }
            }
        }   
    }
}


void  CDisplaySeqalign::x_SetFeatureInfo(CRef<SAlnFeatureInfo> feat_info, 
                                         const CSeq_loc& seqloc, int aln_from, 
                                         int aln_to, int aln_stop, 
                                         char pattern_char, string pattern_id,
                                         string& alternative_feat_str) const
{
    CRef<FeatureInfo> feat(new FeatureInfo);
    feat->seqloc = &seqloc;
    feat->feature_char = pattern_char;
    feat->feature_id = pattern_id;
    
    if(alternative_feat_str != NcbiEmptyString){
        feat_info->feature_string = alternative_feat_str;
    } else {
        //fill feature string
        string line(aln_stop+1, ' ');
        for (int j = aln_from; j <= aln_to; j++){
            line[j] = feat->feature_char;
        }
        feat_info->feature_string = line;
    }
    
    feat_info->aln_range.Set(aln_from, aln_to); 
    feat_info->feature = feat;
}

///add a "|" to the current insert for insert on next rows and return the
///insert end position.
///@param seq: the seq string
///@param insert_aln_pos: the position of insert
///@param aln_start: alnment start position
///@return: the insert end position
///
static int x_AddBar(string& seq, int insert_alnpos, int aln_start){
    int end = (int)seq.size() -1 ;
    int barPos = insert_alnpos - aln_start + 1;
    string addOn;
    if(barPos - end > 1){
        string spacer(barPos - end - 1, ' ');
        addOn += spacer + "|";
    } else if (barPos - end == 1){
        addOn += "|";
    }
    seq += addOn;
    return max<int>((barPos - end), 0);
}


///Add new insert seq to the current insert seq and return the end position of
///the latest insert
///@param cur_insert: the current insert string
///@param new_insert: the new insert string
///@param insert_alnpos: insert position
///@param aln_start: alnment start
///@return: the updated insert end position
///
static int s_AdjustInsert(string& cur_insert, string& new_insert, 
                          int insert_alnpos, int aln_start)
{
    int insertEnd = 0;
    int curInsertSize = (int)cur_insert.size();
    int insertLeftSpace = insert_alnpos - aln_start - curInsertSize + 2;  
    //plus2 because insert is put after the position
    if(curInsertSize > 0){
        _ASSERT(insertLeftSpace >= 2);
    }
    int newInsertSize = (int)new_insert.size();  
    if(insertLeftSpace - newInsertSize >= 1){ 
        //can insert with the end position right below the bar
        string spacer(insertLeftSpace - newInsertSize, ' ');
        cur_insert += spacer + new_insert;
        
    } else { //Need to insert beyond the insert postion
        if(curInsertSize > 0){
            cur_insert += " " + new_insert;
        } else {  //can insert right at the firt position
            cur_insert += new_insert;
        }
    }
    insertEnd = aln_start + (int)cur_insert.size() -1 ; //-1 back to string position
    return insertEnd;
}


void CDisplaySeqalign::x_DoFills(int row, CAlnMap::TSignedRange& aln_range, 
                                 int  aln_start, 
                                 TSInsertInformationList& insert_list, 
                                 list<string>& inserts) const {
    if(!insert_list.empty()){
        string bar(aln_range.GetLength(), ' ');
        
        string seq;
        TSInsertInformationList leftOverInsertList;
        bool isFirstInsert = true;
        int curInsertAlnStart = 0;
        int prvsInsertAlnEnd = 0;
        
        //go through each insert and fills the seq if it can 
        //be filled on the same line.  If not, go to the next line
        NON_CONST_ITERATE(TSInsertInformationList, iter, insert_list) {
            curInsertAlnStart = (*iter)->aln_start;
            //always fill the first insert.  Also fill if there is enough space
            if(isFirstInsert || curInsertAlnStart - prvsInsertAlnEnd >= 1){
                bar[curInsertAlnStart-aln_start+1] = '|';  
                int seqStart = (*iter)->seq_start;
                int seqEnd = seqStart + (*iter)->insert_len - 1;
                string newInsert;
                newInsert = m_AV->GetSeqString(newInsert, row, seqStart,
                                               seqEnd);
                prvsInsertAlnEnd = s_AdjustInsert(seq, newInsert,
                                                  curInsertAlnStart, aln_start);
                isFirstInsert = false;
            } else { //if no space, save the chunk and go to next line 
                bar[curInsertAlnStart-aln_start+1] = '|'; 
                //indicate insert goes to the next line
                prvsInsertAlnEnd += x_AddBar(seq, curInsertAlnStart, aln_start); 
                //May need to add a bar after the current insert sequence 
                //to indicate insert goes to the next line.
                leftOverInsertList.push_back(*iter);    
            }
        }
        //save current insert.  Note that each insert has a bar and sequence
        //below it
        inserts.push_back(bar);
        inserts.push_back(seq);
        //here recursively fill the chunk that don't have enough space
        x_DoFills(row, aln_range, aln_start, leftOverInsertList, inserts);
    }
    
}


void CDisplaySeqalign::x_FillInserts(int row, CAlnMap::TSignedRange& aln_range,
                                     int aln_start, list<string>& inserts,
                                     string& insert_pos_string, 
                                     TSInsertInformationList& insert_list) const
{
    
    string line(aln_range.GetLength(), ' ');
    
    ITERATE(TSInsertInformationList, iter, insert_list){
        int from = (*iter)->aln_start;
        line[from - aln_start + 1] = '\\';
    }
    insert_pos_string = line; 
    //this is the line with "\" right after each insert position
    
    //here fills the insert sequence
    x_DoFills(row, aln_range, aln_start, insert_list, inserts);
}


void CDisplaySeqalign::x_GetInserts(TSInsertInformationList& insert_list,
                                    CAlnMap::TSeqPosList& insert_aln_start, 
                                    CAlnMap::TSeqPosList& insert_seq_start, 
                                    CAlnMap::TSeqPosList& insert_length, 
                                    int line_aln_stop)
{

    while(!insert_aln_start.empty() 
          && (int)insert_aln_start.front() < line_aln_stop){
        CRef<SInsertInformation> insert(new SInsertInformation);
        insert->aln_start = insert_aln_start.front() - 1; 
        //Need to minus one as we are inserting after this position
        insert->seq_start = insert_seq_start.front();
        insert->insert_len = insert_length.front();
        insert_list.push_back(insert);
        insert_aln_start.pop_front();
        insert_seq_start.pop_front();
        insert_length.pop_front();
    }
    
}


string CDisplaySeqalign::x_GetSegs(int row) const 
{
    string segs = NcbiEmptyString;
    if(m_AlignOption & eMergeAlign){ //only show this hsp
        segs = NStr::IntToString(m_AV->GetSeqStart(row))
            + "-" + NStr::IntToString(m_AV->GetSeqStop(row));
    } else { //for all segs
        string idString = m_AV->GetSeqId(1).GetSeqIdString();        
		map<string, struct SAlnLinksParams>::const_iterator iter = m_AlnLinksParams.find(idString);		
		if ( iter != m_AlnLinksParams.end() ){
            segs = iter->second.segs;
        }		
    }
    return segs;
}

	
        
string CDisplaySeqalign::x_GetDumpgnlLink(const list<CRef<CSeq_id> >& ids) const
{
	string dowloadUrl;
    string segs = x_GetSegs(1); //row=1   	
    string label =  CAlignFormatUtil::GetLabel(FindBestChoice(ids, CSeq_id::WorstRank));	
    string url_with_parameters = CAlignFormatUtil::BuildUserUrl(ids, 0, kDownloadUrl,
                                                         m_DbName,
                                                         m_IsDbNa, m_Rid, m_QueryNumber,
                                                         true);
    if (url_with_parameters != NcbiEmptyString) {
        dowloadUrl = CAlignFormatUtil::MapTemplate(kDownloadLink,"download_url",url_with_parameters);
		dowloadUrl = CAlignFormatUtil::MapTemplate(dowloadUrl,"segs",segs);
		dowloadUrl = CAlignFormatUtil::MapTemplate(dowloadUrl,"lnk_displ",kDownloadImg);
		dowloadUrl = CAlignFormatUtil::MapTemplate(dowloadUrl,"label",label);		
    }
    return dowloadUrl;
}


CDisplaySeqalign::DbType CDisplaySeqalign::x_GetDbType(const CSeq_align_set& actual_aln_list) 
{
    //determine if the database has gi by looking at the 1st hit.  
    //Could be wrong but simple for now
    DbType type = eDbTypeNotSet;
    CRef<CSeq_align> first_aln = actual_aln_list.Get().front();
    const CSeq_id& subject_id = first_aln->GetSeq_id(1);
    const CBioseq_Handle& handleTemp  = m_Scope.GetBioseqHandle(subject_id);
    if(handleTemp){
        int giTemp = FindGi(handleTemp.GetBioseqCore()->GetId());
        if (giTemp >0) { 
            type = eDbGi;
        } else if (subject_id.Which() == CSeq_id::e_General){
            const CDbtag& dtg = subject_id.GetGeneral();
            const string& dbName = dtg.GetDb();
            if(NStr::CompareNocase(dbName, "TI") == 0){
                type = eDbGeneral;
            }
        }   
    }
    return type;
}


CRef<CSeq_align_set> 
CDisplaySeqalign::PrepareBlastUngappedSeqalign(const CSeq_align_set& alnset) 
{
    CRef<CSeq_align_set> alnSetRef(new CSeq_align_set);

    ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()){
        const CSeq_align::TSegs& seg = (*iter)->GetSegs();
        if(seg.Which() == CSeq_align::C_Segs::e_Std){
            if(seg.GetStd().size() > 1){ 
                //has more than one stdseg. Need to seperate as each 
                //is a distinct HSP
                ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
                    CRef<CSeq_align> aln(new CSeq_align);
                    if((*iterStdseg)->IsSetScores()){
                        aln->SetScore() = (*iterStdseg)->GetScores();
                    }
                    aln->SetSegs().SetStd().push_back(*iterStdseg);
                    alnSetRef->Set().push_back(aln);
                }
                
            } else {
                alnSetRef->Set().push_back(*iter);
            }
        } else if(seg.Which() == CSeq_align::C_Segs::e_Dendiag){
            if(seg.GetDendiag().size() > 1){ 
                //has more than one dendiag. Need to seperate as each is
                //a distinct HSP
                ITERATE (CSeq_align::C_Segs::TDendiag, iterDendiag,
                         seg.GetDendiag()){
                    CRef<CSeq_align> aln(new CSeq_align);
                    if((*iterDendiag)->IsSetScores()){
                        aln->SetScore() = (*iterDendiag)->GetScores();
                    }
                    aln->SetSegs().SetDendiag().push_back(*iterDendiag);
                    alnSetRef->Set().push_back(aln);
                }
                
            } else {
                alnSetRef->Set().push_back(*iter);
            }
        } else { //Denseg, doing nothing.
            
            alnSetRef->Set().push_back(*iter);
        }
    }
    
    return alnSetRef;
}


bool CDisplaySeqalign::x_IsGeneInfoAvailable(SAlnInfo* aln_vec_info)
{
    const CBioseq_Handle& bsp_handle =
        aln_vec_info->alnvec->GetBioseqHandle(1);
    if (bsp_handle &&
        (m_AlignOption&eHtml) &&
        (m_AlignOption&eLinkout) &&
        (m_AlignOption&eShowGeneInfo))
    {
        CNcbiEnvironment env;
        if (env.Get(GENE_INFO_PATH_ENV_VARIABLE) == kEmptyStr)
        {
            return false;
        }

        const CRef<CBlast_def_line_set> bdlRef 
            =  CAlignFormatUtil::GetBlastDefline(bsp_handle);
        const CBlast_def_line_set::Tdata& bdl = bdlRef->Get();

        ITERATE(CBlast_def_line_set::Tdata, iter, bdl)
        {
            int linkout = m_UseLinkoutDB
                ? CLinkoutDB::GetInstance().GetLinkout(*(*iter)->GetSeqid().front())
                : CAlignFormatUtil::GetLinkout((**iter));
            if (linkout & eGene)
            {
                return true;
            }
        }
    }
    return false;
}


string CDisplaySeqalign::x_GetGeneLinkUrl(int gene_id)
{
    string strGeneLinkUrl = CAlignFormatUtil::GetURLFromRegistry("GENE_INFO");
    AutoPtr<char, ArrayDeleter<char> > buf
        (new char[strGeneLinkUrl.size() + 1024]);
    sprintf(buf.get(), strGeneLinkUrl.c_str(), 
                 gene_id,
                 m_Rid.c_str(),
                 m_IsDbNa ? "nucl" : "prot",
                 m_cur_align);
    strGeneLinkUrl.assign(buf.get());
    return strGeneLinkUrl;
}



string CDisplaySeqalign::x_DisplayGeneInfo(const CBioseq_Handle& bsp_handle,SAlnInfo* aln_vec_info)
{
    CNcbiOstrstream out;
    try
    {
        if (x_IsGeneInfoAvailable(aln_vec_info))
        {
            if (m_GeneInfoReader.get() == 0)
            {
                m_GeneInfoReader.reset(new CGeneInfoFileReader(false));
            }

            int giForGeneLookup = FindGi(bsp_handle.GetBioseqCore()->GetId());

            CGeneInfoFileReader::TGeneInfoList infoList;
            m_GeneInfoReader->GetGeneInfoForGi(giForGeneLookup,infoList);

            CGeneInfoFileReader::TGeneInfoList::const_iterator
                        itInfo = infoList.begin();
            if (itInfo != infoList.end())
                out << "\n";
            for (; itInfo != infoList.end(); itInfo++)
            {
                CRef<CGeneInfo> info = *itInfo;
                string strUrl = x_GetGeneLinkUrl(info->GetGeneId());
                string strInfo;
                info->ToString(strInfo, true, strUrl);
                out << strInfo << "\n";
            }            
        }
    }
    catch (CException& e)
    {
        out << "(Gene info extraction error: "
        << e.GetMsg() << ")" << "\n";
    }
    catch (...)
    {
        out << "(Gene info extraction error)" << "\n";
    }
    string formattedString = CNcbiOstrstreamToString(out);    
    return formattedString;
}

void CDisplaySeqalign::x_DisplayAlignSortInfo(CNcbiOstream& out,string id_label)
{
    string query_buf; 
    map< string, string> parameters_to_change;
    parameters_to_change.insert(map<string, string>::value_type("HSP_SORT", ""));
    CAlignFormatUtil::BuildFormatQueryString(*m_Ctx,parameters_to_change,query_buf);
    out << "\n";
    CAlignFormatUtil::AddSpace(out, 57); 
    out << "Sort alignments for this subject sequence by:\n";
    CAlignFormatUtil::AddSpace(out, 59); 
            
    string hsp_sort_value = m_Ctx->GetRequestValue("HSP_SORT").GetValue();
    int hsp_sort = hsp_sort_value == NcbiEmptyString ? 0 : NStr::StringToInt(hsp_sort_value);
           
    if (hsp_sort != CAlignFormatUtil::eEvalue) {
        out << "<a href=\"Blast.cgi?CMD=Get&" << query_buf 
            << "&HSP_SORT="
            << CAlignFormatUtil::eEvalue
            << "#" << id_label << "\">";
     }
            
     out << "E value";
     if (hsp_sort != CAlignFormatUtil::eEvalue) {
        out << "</a>"; 
     }
           
     CAlignFormatUtil::AddSpace(out, 2);

     if (hsp_sort != CAlignFormatUtil::eScore) {
        out << "<a href=\"Blast.cgi?CMD=Get&" << query_buf 
            << "&HSP_SORT="
            << CAlignFormatUtil::eScore
            << "#" << id_label << "\">";
     }
            
     out << "Score";
     if (hsp_sort != CAlignFormatUtil::eScore) {
        out << "</a>"; 
     }
           
     CAlignFormatUtil::AddSpace(out, 2);

     if (hsp_sort != CAlignFormatUtil::eHspPercentIdentity) {
        out << "<a href=\"Blast.cgi?CMD=Get&" << query_buf 
            << "&HSP_SORT="
            << CAlignFormatUtil::eHspPercentIdentity
            << "#" << id_label << "\">";
     }
     out  << "Percent identity"; 
     if (hsp_sort != CAlignFormatUtil::eHspPercentIdentity) {
        out << "</a>"; 
     }
     out << "\n";
     CAlignFormatUtil::AddSpace(out, 59); 
     if (hsp_sort != CAlignFormatUtil::eQueryStart) {
        out << "<a href=\"Blast.cgi?CMD=Get&" << query_buf 
            << "&HSP_SORT="
            << CAlignFormatUtil::eQueryStart
            << "#" << id_label << "\">";
     } 
     out << "Query start position";
     if (hsp_sort != CAlignFormatUtil::eQueryStart) {
        out << "</a>"; 
     }
     CAlignFormatUtil::AddSpace(out, 2);
           
     if (hsp_sort != CAlignFormatUtil::eSubjectStart) {
        out << "<a href=\"Blast.cgi?CMD=Get&" << query_buf 
            << "&HSP_SORT="
            << CAlignFormatUtil::eSubjectStart
            << "#" << id_label << "\">";
     } 
     out << "Subject start position";
     if (hsp_sort != CAlignFormatUtil::eSubjectStart) {
        out << "</a>"; 
    }
            
    out << "\n";
}

string CDisplaySeqalign::x_FormatAlignSortInfo(string id_label)
{
    string alignSort = m_AlignTemplates->sortInfoTmpl;
    alignSort = CAlignFormatUtil::MapTemplate(alignSort,"id_label",id_label);
    string hsp_sort_value = m_Ctx->GetRequestValue("HSP_SORT").GetValue();
    int hsp_sort = hsp_sort_value == NcbiEmptyString ?  0 : NStr::StringToInt(hsp_sort_value);
    for(size_t  i = 0; i < 5; i++) {
        if(hsp_sort == i) {
            alignSort = CAlignFormatUtil::MapTemplate(alignSort,"sorted_" + NStr::IntToString(hsp_sort),"sortAlnArrowLinkW");                    
        }
        else {
            alignSort = CAlignFormatUtil::MapTemplate(alignSort,"sorted_" + NStr::IntToString(i),""); 
        }
     }            
     return alignSort;
}

void CDisplaySeqalign::x_DisplayBl2SeqLink(CNcbiOstream& out)
{
    const CBioseq_Handle& query_handle=m_AV->GetBioseqHandle(0);
    const CBioseq_Handle& subject_handle=m_AV->GetBioseqHandle(1);
    CSeq_id_Handle query_seqid = GetId(query_handle, eGetId_Best);
    CSeq_id_Handle subject_seqid = GetId(subject_handle, eGetId_Best);
    int query_gi = FindGi(query_handle.GetBioseqCore()->GetId());   
    int subject_gi = FindGi(subject_handle.GetBioseqCore()->GetId());
    
    char buffer[512];
    
    sprintf(buffer, kBl2seqUrl.c_str(), m_Rid.c_str(), 
            query_gi > 0 ? 
            NStr::IntToString(query_gi).c_str():query_seqid.    \
            GetSeqId()->AsFastaString().c_str(),
            subject_gi > 0 ? 
            NStr::IntToString(subject_gi).c_str():subject_seqid. \
            GetSeqId()->AsFastaString().c_str()); 
    out << buffer << "\n";
}

void CDisplaySeqalign::x_DisplayMpvAnchor(CNcbiOstream& out,SAlnInfo* aln_vec_info)
{
    //add id anchor for mapviewer link
    string type_temp = m_BlastType;
    type_temp = NStr::TruncateSpaces(NStr::ToLower(type_temp));
    if(m_AlignOption&eHtml && 
           (type_temp.find("genome") != string::npos ||
            type_temp == "mapview" || 
            type_temp == "mapview_prev" || 
            type_temp == "gsfasta" || type_temp == "gsfasta_prev")){
        string subj_id_str;
        char buffer[126];
        int master_start = m_AV->GetSeqStart(0) + 1;
        int master_stop = m_AV->GetSeqStop(0) + 1;
        int subject_start = m_AV->GetSeqStart(1) + 1;
        int subject_stop = m_AV->GetSeqStop(1) + 1;
    
        m_AV->GetSeqId(1).GetLabel(&subj_id_str, CSeq_id::eContent);
    
        sprintf(buffer, "<a name = %s_%d_%d_%d_%d_%d></a>",
            subj_id_str.c_str(), aln_vec_info->score,
            min(master_start, master_stop),
            max(master_start, master_stop),
            min(subject_start, subject_stop),
            max(subject_start, subject_stop));
    
        out << buffer << "\n"; 
    }
}

string CDisplaySeqalign::x_FormatAlnBlastInfo(SAlnInfo* aln_vec_info)
{
    string evalue_buf, bit_score_buf, total_bit_buf, raw_score_buf;
    CAlignFormatUtil::GetScoreString(aln_vec_info->evalue, 
                                     aln_vec_info->bits, 0, 0, evalue_buf, 
                                     bit_score_buf, total_bit_buf, raw_score_buf);

    string alignParams = m_AlignTemplates->alignInfoTmpl;
    if (m_SeqalignSetRef->Get().front()->CanGetType() && 
           m_SeqalignSetRef->Get().front()->GetType() == CSeq_align_Base::eType_global)
    {
        //out<<" NW Score = "<< aln_vec_info->score; ///??? Add NW score               
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_score",aln_vec_info->score);
    }
    else
    {
        //out<<" Score = "<<bit_score_buf<<" ";        
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_score",bit_score_buf);
        //out<<"bits ("<<aln_vec_info->score<<"),"<<"  ";
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_score_bits",aln_vec_info->score);
            
        //out<<"Expect";
        //out << " = " << evalue_buf;
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_eval",evalue_buf);
        if (aln_vec_info->sum_n > 0) {
            //out << "(" << aln_vec_info->sum_n << ")";///???SumN - get rid - check with Tom
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_sumN",aln_vec_info->sum_n);
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"sumNshow","shown");
        }
        else {
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_sumN","");
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"sumNshow","");
        }
                
        if (aln_vec_info->comp_adj_method == 1){
            //out << ", Method: Composition-based stats.";
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth","Composition-based stats.");
            alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth_hide","");//???? is that the same for all aligns??? 
        }
        else if (aln_vec_info->comp_adj_method == 2){
           //out << ", Method: Compositional matrix adjust.";
           alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth","Compositional matrix adjust.");
           alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth_hide","");//???? is that the same for all aligns??? 
        }
        else {
          alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth_hide","hidden");//???? is that the same for all aligns??? 
          alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_meth","");
        }        
    }    
    return alignParams;    
}
//sumN - hidden, cbs_md - shown, aln_frame - hidden



void CDisplaySeqalign::x_DisplayAlignInfo(CNcbiOstream& out, 
                                           SAlnInfo* aln_vec_info)
{
    string evalue_buf, bit_score_buf, total_bit_buf, raw_score_buf;
    CAlignFormatUtil::GetScoreString(aln_vec_info->evalue, 
                                         aln_vec_info->bits, 0, 0, evalue_buf, 
                                         bit_score_buf, total_bit_buf, raw_score_buf);
        
    CRef<CSeq_align> first_aln = m_SeqalignSetRef->Get().front();
    if (m_SeqalignSetRef->Get().front()->CanGetType() && 
       m_SeqalignSetRef->Get().front()->GetType() == CSeq_align_Base::eType_global)
    {
        out<<" NW Score = "<< aln_vec_info->score;
    }
    else
    {
        out<<" Score = "<<bit_score_buf<<" ";
        out<<"bits ("<<aln_vec_info->score<<"),"<<"  ";
        out<<"Expect";
        if (aln_vec_info->sum_n > 0) {
        out << "(" << aln_vec_info->sum_n << ")";
        }
        out << " = " << evalue_buf;
        if (aln_vec_info->comp_adj_method == 1)
        out << ", Method: Composition-based stats.";
        else if (aln_vec_info->comp_adj_method == 2)
        out << ", Method: Compositional matrix adjust.";
    }
    out << "\n";
}

//1. Display defline(s)           
//2. Display Gene info
//3. Display Bl2Seq TBLASTX link
//4. add id anchor for mapviewer link
void CDisplaySeqalign::x_ShowAlnvecInfo(CNcbiOstream& out, 
                                           SAlnInfo* aln_vec_info,
                                           bool show_defline) 
{
	bool showSortControls = false;
    if(show_defline) {        
		const CBioseq_Handle& bsp_handle=m_AV->GetBioseqHandle(1); 
		if(m_AlignOption&eShowBlastInfo) {
			if(!(m_AlignOption & eShowNoDeflineInfo)){
				//1. Display defline(s),Gene info				
				string deflines = x_PrintDefLine(bsp_handle, aln_vec_info);
				out<< deflines;
				//2. Format Gene info
				string geneInfo = x_DisplayGeneInfo(bsp_handle,aln_vec_info);
				out<< geneInfo;							       
			}       
        
			if((m_AlignOption&eHtml) && (m_AlignOption&eShowBlastInfo)
				&& (m_AlignOption&eShowBl2seqLink)) {
				//3. Display Bl2Seq TBLASTX link
				x_DisplayBl2SeqLink(out);
			}
			out << "\n";
		}
        showSortControls = true;
    }
    if (m_AlignOption&eShowBlastInfo) {
        //4. add id anchor for mapviewer link
        x_DisplayMpvAnchor(out,aln_vec_info);    
    }
    
    //Displays sorting controls, features, Score, Expect, Idnt,Gaps,strand,positives,frames etc
    x_DisplaySingleAlignParams(out, aln_vec_info,showSortControls);
    string alignRows = x_DisplayRowData(aln_vec_info->alnRowInfo);
    out << alignRows;            
}



//fill one defline info, using  <@ALN_DEFLINE_ROW@>
string
CDisplaySeqalign::x_MapDefLine(SAlnDispParams *alnDispParams,bool isFirst, bool linkout,bool hideDefline)
{
	/*
    string firstSeqClassInfo = (isFirst) ? "" : "hidden"; //hide ">" sign if not first seq align	
	string alnDefLine  = CAlignFormatUtil::MapTemplate(m_AlignTemplates->alnDefLineTmpl,"alnSeqSt",firstSeqClassInfo);
	*/
	string alnDflnShow = (hideDefline) ? "hidden" : "";	
	string alnDefLine = CAlignFormatUtil::MapTemplate(m_AlignTemplates->alnDefLineTmpl,"alnDflnShow",alnDflnShow);

	string alnGi = (m_AlignOption&eShowGi && alnDispParams->gi > 0) ? "gi|" + NStr::IntToString(alnDispParams->gi) + "|" : "";
	string seqid;					
    if(!(alnDispParams->seqID->AsFastaString().find("gnl|BL_ORD_ID") != string::npos)){							 
		seqid = alnDispParams->seqID->AsFastaString();        
    }
	
	if(alnDispParams->id_url != NcbiEmptyString) {
		string seqInfo  = CAlignFormatUtil::MapTemplate(m_AlignTemplates->alnSeqInfoTmpl,"aln_url",alnDispParams->id_url);
		string trgt = (m_AlignOption & eNewTargetWindow) ? "TARGET=\"EntrezView\"" : "";

		seqInfo = CAlignFormatUtil::MapTemplate(seqInfo,"aln_target",trgt);
		alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"seq_info",seqInfo);        		
		alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"aln_gi",alnGi);    
		alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"aln_seqid",seqid);               
    }
	else {
		alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"seq_info",alnGi + seqid); 
	}
	string alnIdLbl = (alnDispParams->gi != 0) ? NStr::IntToString(alnDispParams->gi) : alnDispParams->seqID->GetSeqIdString();
	alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"alnIdLbl",alnIdLbl);
	string linkoutStr, dnldLinkStr;
	if (linkout) {
		linkoutStr = (!alnDispParams->linkoutStr.empty()) ? alnDispParams->linkoutStr : "";		
		dnldLinkStr = alnDispParams->dumpGnlUrl;
	}		
	alnDefLine  = CAlignFormatUtil::MapTemplate(alnDefLine ,"alnLinkout",linkoutStr);
	alnDefLine  = CAlignFormatUtil::MapTemplate(alnDefLine ,"dndlLinkt",dnldLinkStr);	
	alnDefLine = CAlignFormatUtil::MapTemplate(alnDefLine,"alnTitle",alnDispParams->title);			
	return alnDefLine;
}

string
CDisplaySeqalign::x_FormatDefLinesHeader(const CBioseq_Handle& bsp_handle,SAlnInfo* aln_vec_info)
{
    CNcbiOstrstream out;    
    string deflines,accession,id_label;
    if(bsp_handle){        
        const CRef<CBlast_def_line_set> bdlRef =  CAlignFormatUtil::GetBlastDefline(bsp_handle);
        const list< CRef< CBlast_def_line > >& bdl = bdlRef->Get();
        bool isFirst = true;
        int firstGi = 0;
		m_NumBlastDefLines = 0;
        m_cur_align++;		
		SAlnDispParams *alnDispParams;
        if(bdl.empty()){ //no blast defline struct, should be no such case now
            //actually not so fast...as we now fetch from entrez even when it's not in blast db
            //there is no blast defline in such case.
			alnDispParams = x_FillAlnDispParams(bsp_handle);
			string alnDefLine = x_MapDefLine(alnDispParams,isFirst,false,false);
		    id_label = (alnDispParams->gi != 0) ? NStr::IntToString(alnDispParams->gi) : alnDispParams->label;
			accession = alnDispParams->label;
			delete alnDispParams;
			deflines = alnDefLine;
        } else {
            //format each defline 
            int numBdl = 0;
            for(list< CRef< CBlast_def_line > >::const_iterator 
                    iter = bdl.begin(); iter != bdl.end(); iter++){
                numBdl++;
				alnDispParams = x_FillAlnDispParams(*iter,aln_vec_info->use_this_gi,firstGi,bsp_handle.GetBioseqCore()->IsNa(),(int)bsp_handle.GetBioseqLength());
				if(alnDispParams) {
					bool hideDefline = (int)(bdl.size()) > k_MaxDeflinesToShow && numBdl >= k_MinDeflinesToShow + 1;
					string alnDefLine = x_MapDefLine(alnDispParams,isFirst,m_AlignOption&eLinkout,hideDefline);                    
                    if(isFirst){
                        firstGi = alnDispParams->gi;
						id_label = (alnDispParams->gi != 0) ? NStr::IntToString(alnDispParams->gi) : alnDispParams->label;
                        accession = alnDispParams->label;
                    }                    
					deflines += alnDefLine;	
                    isFirst = false;					
					delete alnDispParams;
                }
            }			
			m_NumBlastDefLines = numBdl;            
        }
    }    
	//fill deflines
	string alignInfo = CAlignFormatUtil::MapTemplate(m_AlignTemplates->alignHeaderTmpl,"aln_deflines",deflines);

	//fill length
	int seqLength = bsp_handle.GetBioseqLength();
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnSeqLength", NStr::IntToString(seqLength));

	//fill multiple titles
	int alnSeqTitlesNum = (m_NumBlastDefLines > k_MaxDeflinesToShow) ? m_NumBlastDefLines - k_MinDeflinesToShow : 0;
	string alnSeqTitlesShow = (m_NumBlastDefLines > k_MaxDeflinesToShow) ? "" : "hidden";
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnSeqTitlesNum", NStr::IntToString(alnSeqTitlesNum));
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnSeqTitlesShow",alnSeqTitlesShow);				

	//fill id info
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnID",id_label);				
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"firstSeqID",accession);	
   
	//fill sequence checkbox
	string seqRetrieval = ((m_AlignOption&eSequenceRetrieval) && m_CanRetrieveSeq) ? "" : "hidden";
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnSeqGi",id_label);
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnQueryNum",NStr::IntToString(m_QueryNumber));				
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"alnSeqRet",seqRetrieval);

	//fill sort info
	string sortInfo;	
    if(m_AlnLinksParams[m_AV->GetSeqId(1).GetSeqIdString()].hspNumber > 1 &&	
		m_AlignOption & eShowSortControls){
		//3. Display sort info
		sortInfo = x_FormatAlignSortInfo(id_label);					
	}
	alignInfo  = CAlignFormatUtil::MapTemplate(alignInfo,"sortInfo",sortInfo);
	
    return alignInfo;
}


//1. Display defline(s)           
//2. Display Gene info
//3. Display Bl2Seq TBLASTX link
void CDisplaySeqalign::x_ShowAlnvecInfoTemplate(CNcbiOstream& out, 
                                           SAlnInfo* aln_vec_info,
                                           bool show_defline,
                                           bool showSortControls) 
{
    if(show_defline) {        
        const CBioseq_Handle& bsp_handle=m_AV->GetBioseqHandle(1); 
		if(m_AlignOption&eShowBlastInfo) {
			if(!(m_AlignOption & eShowNoDeflineInfo)){            
				//1. Display defline(s),Gene info
				string alignHeader = x_FormatDefLinesHeader(bsp_handle, aln_vec_info);
				//2. Format Gene info
				string geneInfo = x_DisplayGeneInfo(bsp_handle,aln_vec_info);				
                alignHeader = CAlignFormatUtil::MapTemplate(alignHeader,"aln_gene_info",geneInfo); 
                out<< alignHeader;			            
			}       
        
			if((m_AlignOption&eHtml) && (m_AlignOption&eShowBlastInfo)
				&& (m_AlignOption&eShowBl2seqLink)) {
				//3. Display Bl2Seq TBLASTX link
				x_DisplayBl2SeqLink(out);
			}
			out << "\n";
		}        
		//start counting hsp
		m_currAlignHsp = 0;
    }
    if (m_AlignOption&eShowBlastInfo) {
        //4. add id anchor for mapviewer link
        x_DisplayMpvAnchor(out,aln_vec_info);    
    }
    
	if (m_AlignOption&eShowBlastInfo && showSortControls && 
           m_AlnLinksParams[m_AV->GetSeqId(1).GetSeqIdString()].hspNumber > 1 &&
									m_AlignOption & eShowSortControls){
       //5. Display sort info
       string sortInfo = x_FormatAlignSortInfo(aln_vec_info->id_label);
	   out << sortInfo;    
    }

    //Displays sorting controls, features, Score, Expect, Idnt,Gaps,strand,positives,frames etc    
    string alignInfo = x_FormatSingleAlign(aln_vec_info);
    out << alignInfo;        
}

void CDisplaySeqalign::x_DisplayAlnvecInfo(CNcbiOstream& out, 
                                           SAlnInfo* aln_vec_info,
                                           bool show_defline,
                                           bool showSortControls) 
{

    m_AV = aln_vec_info->alnvec;
	//Calculate Dynamic Features in aln_vec_info    
    x_PrepareDynamicFeatureInfo(aln_vec_info);    
    //Calculate row data for actual alignment display
    aln_vec_info->alnRowInfo = x_PrepareRowData();

    //Calculate indentity data in aln_vec_info  
    if((m_AlignOption & eShowBlastInfo) || (m_AlignOption & eShowMiddleLine)){
        x_PrepareIdentityInfo(aln_vec_info);
    }
	if(!m_AlignTemplates) {
		x_ShowAlnvecInfo(out,aln_vec_info,show_defline);
	}
	else {
		x_ShowAlnvecInfoTemplate(out,aln_vec_info,show_defline,showSortControls);
	}      
    
    delete aln_vec_info->alnRowInfo;

    out<<"\n";
}


//Displays features, Score Expect, Idnt,Gaps,strand
void CDisplaySeqalign::x_DisplaySingleAlignParams(CNcbiOstream& out, 
                                                SAlnInfo* aln_vec_info,
                                                bool showSortControls)  
{
    if (m_AlignOption&eShowBlastInfo) {            
    
        if(showSortControls && m_AlignOption&eHtml &&
           m_AlnLinksParams[m_AV->GetSeqId(1).GetSeqIdString()].hspNumber > 1 &&
           m_AlignOption & eShowSortControls){
            //3. Display sort info
            x_DisplayAlignSortInfo(out,aln_vec_info->id_label);
        }

        //output dynamic feature lines
        if(aln_vec_info->feat_list.size() > 0 || aln_vec_info->feat5 || aln_vec_info->feat3 ){        
            //6. Display Dynamic Features
            x_PrintDynamicFeatures(out,aln_vec_info);             
        }
        
        //7. Display score,bits,expect,method
        x_DisplayAlignInfo(out,aln_vec_info);
    }
    
    if((m_AlignOption & eShowBlastInfo) || (m_AlignOption & eShowMiddleLine)){
        //8.Display Identities,positives,strand, frames etc
        //x_DisplayIdentityInfo(aln_vec_info->alnRowInfo, out);
        s_DisplayIdentityInfo(out, 
                              (int)m_AV->GetAlnStop(), 
                              aln_vec_info->identity, 
                              aln_vec_info->positive, 
                              aln_vec_info->match, 
                              aln_vec_info->gap,
                              m_AV->StrandSign(0), 
                              m_AV->StrandSign(1),
                              aln_vec_info->alnRowInfo->frame[0], 
                              aln_vec_info->alnRowInfo->frame[1], 
                              ((m_AlignType & eProt) != 0 ? true : false));
    }   
}

//Displays features, Score Expect, Idnt,Gaps,strand
string CDisplaySeqalign::x_FormatSingleAlign(SAlnInfo* aln_vec_info)                                                  
{
    string alignInfo;
    if (m_AlignOption&eShowBlastInfo) { 
        
        //7. Display score,bits,expect,method
        alignInfo = x_FormatAlnBlastInfo(aln_vec_info);

        //8.Display Identities,positives,strands, frames etc        
        alignInfo = x_FormatIdentityInfo(alignInfo, aln_vec_info);
        
        //output dynamic feature lines
        //only for aln_vec_info->feat_list.size() > 0 || aln_vec_info->feat5 || aln_vec_info->feat3
        //6. Display Dynamic Features            
        alignInfo = x_FormatDynamicFeaturesInfo(alignInfo, aln_vec_info);    
    }
	m_currAlignHsp++;
	string alignRowsTemplate = (m_currAlignHsp == m_AlnLinksParams[m_AV->GetSeqId(1).GetSeqIdString()].hspNumber) ? m_AlignTemplates->alignRowTmplLast : m_AlignTemplates->alignRowTmpl;
	
    string alignRows = x_DisplayRowData(aln_vec_info->alnRowInfo);
    alignRows = CAlignFormatUtil::MapTemplate(alignRowsTemplate,"align_rows",alignRows);
    alignInfo += alignRows;
    return alignInfo;
    //if((m_AlignOption & eShowBlastInfo) || (m_AlignOption & eShowMiddleLine)){///ignore eShowMiddleLine for templates
    //}   
}




void CDisplaySeqalign::x_PrepareDynamicFeatureInfo(SAlnInfo* aln_vec_info)
{
    aln_vec_info->feat5 = NULL;
    aln_vec_info->feat3 = NULL;
    aln_vec_info->feat_list.clear();    
    //Calculate Dynamic Features in aln_vec_info               
    if((m_AlignOption&eDynamicFeature) 
        && (int)m_AV->GetBioseqHandle(1).GetBioseqLength() 
        >= k_GetDynamicFeatureSeqLength){
        if(m_DynamicFeature){
            const CSeq_id& subject_seqid = m_AV->GetSeqId(1);            
            const CRange<TSeqPos>& range = m_AV->GetSeqRange(1);
            aln_vec_info->actual_range = range;
            if(range.GetFrom() > range.GetTo()){
                aln_vec_info->actual_range.Set(range.GetTo(), range.GetFrom());
            }
            string id_str;
            subject_seqid.GetLabel(&id_str, CSeq_id::eBoth);
            const CBioseq_Handle& subject_handle=m_AV->GetBioseqHandle(1);
            aln_vec_info->subject_gi = FindGi(subject_handle.GetBioseqCore()->GetId());
            aln_vec_info->feat_list  =  m_DynamicFeature->GetFeatInfo(id_str, aln_vec_info->actual_range, aln_vec_info->feat5, aln_vec_info->feat3, 2);
        } 
    }
}

static string s_MapFeatureURL(string viewerURL,                           
                              int subject_gi,  
                              string db,                              
                              int fromRange, 
                              int toRange,
                              string rid)
{    
    string url_link = CAlignFormatUtil::MapTemplate(viewerURL,"db",db);
    url_link = CAlignFormatUtil::MapTemplate(url_link,"gi",subject_gi);    
    url_link = CAlignFormatUtil::MapTemplate(url_link,"rid",rid); 
    url_link = CAlignFormatUtil::MapTemplate(url_link,"from",fromRange); 
    url_link = CAlignFormatUtil::MapTemplate(url_link,"to",toRange); 
    return url_link;
}

string CDisplaySeqalign::x_FormatOneDynamicFeature(string viewerURL,
                                                   int subject_gi,                                                    
                                                   int fromRange, 
                                                   int toRange,
                                                   string featText)
{
    string alignFeature = m_AlignTemplates->alignFeatureTmpl;
    if(subject_gi > 0){                   
        alignFeature = CAlignFormatUtil::MapTemplate(alignFeature,"aln_feat_info",m_AlignTemplates->alignFeatureLinkTmpl);
        string url = s_MapFeatureURL(viewerURL,
                                     subject_gi,
                                     string(m_IsDbNa ? "nucleotide" : "protein"),
                                     fromRange + 1,
                                     toRange + 1,
                                     m_Rid);
        alignFeature = CAlignFormatUtil::MapTemplate(alignFeature,"aln_feat_url",url);
        alignFeature = CAlignFormatUtil::MapTemplate(alignFeature,"aln_feat",featText);
    }
    else {
        alignFeature = CAlignFormatUtil::MapTemplate(alignFeature,"aln_feat_info",featText);
    }
    return alignFeature;
}


//6. Display Dynamic Features
string CDisplaySeqalign::x_FormatDynamicFeaturesInfo(string alignInfo, SAlnInfo* aln_vec_info) 
{
    string alignParams = alignInfo;
    //string alignFeature = m_AlignTemplates->alignFeatureTmpl;
    
    //char urlBuf[2048];   
    string viewerURL = CAlignFormatUtil::GetURLFromRegistry("ENTREZ_SUBSEQ_TM");

    string allAlnFeatures = "";
    if(aln_vec_info->feat_list.size() > 0) { //has feature in this range
        ITERATE(vector<SFeatInfo*>, iter, aln_vec_info->feat_list){            
                
            string alignFeature = x_FormatOneDynamicFeature(viewerURL,
                                                     aln_vec_info->subject_gi,
                                                     (*iter)->range.GetFrom(),
                                                     (*iter)->range.GetTo(),
                                                     (*iter)->feat_str);
                                                   
            ///TO DO: NO hyperlink if aln_vec_info->subject_gi == 0            
            
            allAlnFeatures += alignFeature;            
        }
    } else {  //show flank features
        if(aln_vec_info->feat5 || aln_vec_info->feat3){   
            //TO DO: Check if we need that
            //out << " Features flanking this part of subject sequence:" << "\n";
        }
        if(aln_vec_info->feat5){            
            string alignFeature = x_FormatOneDynamicFeature(viewerURL,
                                                     aln_vec_info->subject_gi,
                                                     aln_vec_info->feat5->range.GetFrom(),
                                                     aln_vec_info->feat5->range.GetTo(),
                                                     NStr::IntToString(aln_vec_info->actual_range.GetFrom() - aln_vec_info->feat5->range.GetTo()) + (string)" bp at 5' side: " + aln_vec_info->feat5->feat_str);
            allAlnFeatures += alignFeature;
        }
        if(aln_vec_info->feat3){     
            
            string alignFeature = x_FormatOneDynamicFeature(viewerURL,
                                                     aln_vec_info->subject_gi,
                                                     aln_vec_info->feat3->range.GetFrom(),
                                                     aln_vec_info->feat3->range.GetTo(),
                                                     NStr::IntToString(aln_vec_info->feat3->range.GetFrom() - aln_vec_info->actual_range.GetTo()) + (string)" bp at 3' side: " + aln_vec_info->feat3->feat_str);
            allAlnFeatures += alignFeature;
        }
    }
    if(!allAlnFeatures.empty()) {        
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"all_aln_features",allAlnFeatures);
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_feat_show","");        
    }
    else {
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"all_aln_features","");
        alignParams = CAlignFormatUtil::MapTemplate(alignParams,"aln_feat_show","hidden");
    }  
    return alignParams;
}

void CDisplaySeqalign::x_PrintDynamicFeatures(CNcbiOstream& out,SAlnInfo* aln_vec_info) 
{
    string l_EntrezSubseqUrl = CAlignFormatUtil::GetURLFromRegistry("ENTREZ_SUBSEQ");

    if(aln_vec_info->feat_list.size() > 0) { //has feature in this range
        out << " Features in this part of subject sequence:" << "\n";
        ITERATE(vector<SFeatInfo*>, iter, aln_vec_info->feat_list){
            out << "   ";
            if(m_AlignOption&eHtml && aln_vec_info->subject_gi > 0){
                string featStr = s_MapFeatureURL(l_EntrezSubseqUrl, 
                                              aln_vec_info->subject_gi,
                                              m_IsDbNa ? "nucleotide" : "protein",  
                                              (*iter)->range.GetFrom() +1 , 
                                              (*iter)->range.GetTo() + 1,
                                              m_Rid);                
                out << featStr;
            }  
            out << (*iter)->feat_str;
            if(m_AlignOption&eHtml && aln_vec_info->subject_gi > 0){
                out << "</a>";
            }  
            out << "\n";
        }
    } else {  //show flank features
        if(aln_vec_info->feat5 || aln_vec_info->feat3){   
            out << " Features flanking this part of subject sequence:" << "\n";
        }
        if(aln_vec_info->feat5){
            out << "   ";
            if(m_AlignOption&eHtml && aln_vec_info->subject_gi > 0){
                string featStr = s_MapFeatureURL(l_EntrezSubseqUrl, 
                                              aln_vec_info->subject_gi,
                                              m_IsDbNa ? "nucleotide" : "protein",  
                                              aln_vec_info->feat5->range.GetFrom() + 1 , 
                                              aln_vec_info->feat5->range.GetTo() + 1,
                                              m_Rid);
                out << featStr;
            }  
            out << aln_vec_info->actual_range.GetFrom() - aln_vec_info->feat5->range.GetTo() 
                << " bp at 5' side: " << aln_vec_info->feat5->feat_str;
            if(m_AlignOption&eHtml && aln_vec_info->subject_gi > 0){
                out << "</a>";
            }  
            out << "\n";
        }
        if(aln_vec_info->feat3){
            out << "   ";
            if(m_AlignOption&eHtml && aln_vec_info->subject_gi > 0){
                string featStr = s_MapFeatureURL(l_EntrezSubseqUrl, 
                                              aln_vec_info->subject_gi,
                                              m_IsDbNa ? "nucleotide" : "protein",  
                                              aln_vec_info->feat3->range.GetFrom() + 1 , 
                                              aln_vec_info->feat3->range.GetTo() + 1,
                                              m_Rid);
                out << featStr;
            }
            out << aln_vec_info->feat3->range.GetFrom() - aln_vec_info->actual_range.GetTo() 
                << " bp at 3' side: " << aln_vec_info->feat3->feat_str;
            if(m_AlignOption&eHtml){
                out << "</a>";
            }  
            out << "\n";
        }
    }
    if(aln_vec_info->feat_list.size() > 0 || aln_vec_info->feat5 || aln_vec_info->feat3 ){
        out << "\n";
    }
}

void 
CDisplaySeqalign::x_FillLocList(TSAlnSeqlocInfoList& loc_list, 
                        const list< CRef<CSeqLocInfo> >* masks) const
{
    if ( !masks ) {
        return;
    }

    ITERATE(TMaskedQueryRegions, iter, *masks) {
        CRef<SAlnSeqlocInfo> alnloc(new SAlnSeqlocInfo);
        bool has_valid_loc = false;
        for (int i=0; i<m_AV->GetNumRows(); i++){
            TSeqRange loc_range(**iter);
            if((*iter)->GetInterval().GetId().Match(m_AV->GetSeqId(i)) &&
               m_AV->GetSeqRange(i).IntersectingWith(loc_range)){
                int actualAlnStart = 0, actualAlnStop = 0;
                if(m_AV->IsPositiveStrand(i)){
                    actualAlnStart =
                        m_AV->GetAlnPosFromSeqPos(i, 
                                                  (*iter)->GetInterval().GetFrom(),
                                                          CAlnMap::eBackwards, true);
                    actualAlnStop =
                        m_AV->GetAlnPosFromSeqPos(i, 
                                                  (*iter)->GetInterval().GetTo(),
                                                  CAlnMap::eBackwards, true);
                } else {
                    actualAlnStart =
                        m_AV->GetAlnPosFromSeqPos(i, 
                                                  (*iter)->GetInterval().GetTo(),
                                                  CAlnMap::eBackwards, true);
                    actualAlnStop =
                        m_AV->GetAlnPosFromSeqPos(i, 
                                                  (*iter)->GetInterval().GetFrom(),
                                                  CAlnMap::eBackwards, true);
                }
                alnloc->aln_range.Set(actualAlnStart, actualAlnStop);  
                has_valid_loc = true;
                break;
            }
        }
        if (has_valid_loc) {
            alnloc->seqloc = *iter;   
            loc_list.push_back(alnloc);
        }
    }
}


void
CDisplaySeqalign::x_GetQueryFeatureList(int row_num, int aln_stop,
                                        vector<TSAlnFeatureInfoList>& retval) 
                                        const
{
    retval.clear();
    retval.resize(row_num);
    //list<SAlnFeatureInfo*>* bioseqFeature= new list<SAlnFeatureInfo*>[row_num];
    if(m_QueryFeature){
        for (list<FeatureInfo*>::iterator iter=m_QueryFeature->begin(); 
             iter!=m_QueryFeature->end(); iter++){
            for(int i = 0; i < row_num; i++){
                if((*iter)->seqloc->GetInt().GetId().Match(m_AV->GetSeqId(i))){
                    int actualSeqStart = 0, actualSeqStop = 0;
                    if(m_AV->IsPositiveStrand(i)){
                        if((*iter)->seqloc->GetInt().GetFrom() 
                           < m_AV->GetSeqStart(i)){
                            actualSeqStart = m_AV->GetSeqStart(i);
                        } else {
                            actualSeqStart = (*iter)->seqloc->GetInt().GetFrom();
                        }
                        
                        if((*iter)->seqloc->GetInt().GetTo() >
                           m_AV->GetSeqStop(i)){
                            actualSeqStop = m_AV->GetSeqStop(i);
                        } else {
                            actualSeqStop = (*iter)->seqloc->GetInt().GetTo();
                        }
                    } else {
                        if((*iter)->seqloc->GetInt().GetFrom() 
                           < m_AV->GetSeqStart(i)){
                            actualSeqStart = (*iter)->seqloc->GetInt().GetFrom();
                        } else {
                            actualSeqStart = m_AV->GetSeqStart(i);
                        }
                        
                        if((*iter)->seqloc->GetInt().GetTo() > 
                           m_AV->GetSeqStop(i)){
                            actualSeqStop = (*iter)->seqloc->GetInt().GetTo();
                        } else {
                            actualSeqStop = m_AV->GetSeqStop(i);
                        }
                    }
                    int alnFrom = m_AV->GetAlnPosFromSeqPos(i, actualSeqStart);
                    int alnTo = m_AV->GetAlnPosFromSeqPos(i, actualSeqStop);
                    
                    CRef<SAlnFeatureInfo> featInfo(new SAlnFeatureInfo);
                    string tempFeat = NcbiEmptyString;
                    x_SetFeatureInfo(featInfo, *((*iter)->seqloc), alnFrom, 
                                     alnTo,  aln_stop, (*iter)->feature_char,
                                     (*iter)->feature_id, tempFeat);    
                    retval[i].push_back(featInfo);
                }
            }
        }
    }
}


void CDisplaySeqalign::x_FillSeqid(string& id, int row) const
{
    static string kQuery("Query");
    static string kSubject("Sbjct");

#ifdef CTOOLKIT_COMPATIBLE
    /* Facilitates comparing formatted output using diff */
    static bool value_set = false;
    if ( !value_set ) {
        if (getenv("CTOOLKIT_COMPATIBLE")) {
            kQuery.append(":");
            kSubject.append(":");
        }
        value_set = true;
    }
#endif /* CTOOLKIT_COMPATIBLE */

    if(m_AlignOption & eShowBlastStyleId) {
        if(row==0){//query
            id=kQuery;
        } else {//hits
            if (!(m_AlignOption&eMergeAlign)){
                //hits for pairwise 
                id=kSubject;
            } else {
                if(m_AlignOption&eShowGi){
                    int gi = 0;
                    if(m_AV->GetSeqId(row).Which() == CSeq_id::e_Gi){
                        gi = m_AV->GetSeqId(row).GetGi();
                    }
                    if(!(gi > 0)){
                        gi = s_GetGiForSeqIdList(m_AV->GetBioseqHandle(row).\
                                                 GetBioseqCore()->GetId());
                    }
                    if(gi > 0){
                        id=NStr::IntToString(gi);
                    } else {
                        const CRef<CSeq_id> wid 
                            = FindBestChoice(m_AV->GetBioseqHandle(row).\
                                             GetBioseqCore()->GetId(), 
                                             CSeq_id::WorstRank);
                        id = CAlignFormatUtil::GetLabel(wid).c_str();
                    }
                } else {
                    const CRef<CSeq_id> wid 
                        = FindBestChoice(m_AV->GetBioseqHandle(row).\
                                         GetBioseqCore()->GetId(), 
                                         CSeq_id::WorstRank);
                    id = CAlignFormatUtil::GetLabel(wid).c_str();
                }           
            }
        }
    } else {
        if(m_AlignOption&eShowGi){
            int gi = 0;
            if(m_AV->GetSeqId(row).Which() == CSeq_id::e_Gi){
                gi = m_AV->GetSeqId(row).GetGi();
            }
            if(!(gi > 0)){
                gi = s_GetGiForSeqIdList(m_AV->GetBioseqHandle(row).\
                                         GetBioseqCore()->GetId());
            }
            if(gi > 0){
                id=NStr::IntToString(gi);
            } else {
                const CRef<CSeq_id> wid 
                    = FindBestChoice(m_AV->GetBioseqHandle(row).\
                                     GetBioseqCore()->GetId(),
                                     CSeq_id::WorstRank);
                id = CAlignFormatUtil::GetLabel(wid).c_str();
            }
        } else {
            const CRef<CSeq_id> wid 
                = FindBestChoice(m_AV->GetBioseqHandle(row).\
                                 GetBioseqCore()->GetId(), 
                                 CSeq_id::WorstRank);
            id = CAlignFormatUtil::GetLabel(wid).c_str();
        }     
    }
}



void CDisplaySeqalign::x_PreProcessSeqAlign(CSeq_align_set &actual_aln_list)
{
    int num_align = 0;
    //get segs first and get hspNumber,segs and subjRange per sequence in alignment
	string toolUrl = NcbiEmptyString;
	if(m_AlignOption & eHtml){
		toolUrl = m_Reg->Get(m_BlastType, "TOOL_URL");
	}
	if(toolUrl.find("report=graph") != string::npos || //do it for query anchored also for report=graph
       (!(m_AlignOption & eMergeAlign) &&
	     (toolUrl.find("dumpgnl.cgi") != string::npos			
		  || (m_AlignOption & eLinkout)
		  || (m_AlignOption & eHtml && m_AlignOption & eShowBlastInfo)))) {
			/*need to construct segs for dumpgnl and
			get sub-sequence for long sequences*/
		
        for (CSeq_align_set::Tdata::const_iterator
                     iter =  actual_aln_list.Get().begin();
                 iter != actual_aln_list.Get().end()
                     && num_align<m_NumAlignToShow; iter++, num_align++) {

			CConstRef<CSeq_id> subid;
            subid = &((*iter)->GetSeq_id(1));
            string idString = subid->GetSeqIdString();
			
			x_CalcUrlLinksParams(**iter,idString,toolUrl);//sets m_AlnLinksParams->segs,hspNum, subjRange 
        }
    }
}


void CDisplaySeqalign::x_CalcUrlLinksParams(const CSeq_align& align, string idString,string toolUrl)
{
    //make alnvector	
    CRef<CAlnVec> avRef = x_GetAlnVecForSeqalign(align);    
	
	bool first = m_AlnLinksParams.count(idString) == 0;    
    struct SAlnLinksParams *alnLinksParam = first ? new SAlnLinksParams : &m_AlnLinksParams[idString];
	

	if (toolUrl.find("dumpgnl.cgi") != string::npos || (m_AlignOption & eLinkout))  {
		if(!first){
			alnLinksParam->segs += ",";
		}
		alnLinksParam->segs += NStr::IntToString(avRef->GetSeqStart(1))
							   + "-" +
							   NStr::IntToString(avRef->GetSeqStop(1));        
	}
    
	if (toolUrl.find("report=graph") != string::npos) {
		TSeqPos from = (avRef->GetSeqStart(1)> avRef->GetSeqStop(1)) ? avRef->GetSeqStop(1) : avRef->GetSeqStart(1);
        TSeqPos to =   (avRef->GetSeqStart(1)> avRef->GetSeqStop(1)) ? avRef->GetSeqStart(1) : avRef->GetSeqStop(1);			
		if(first) {
			alnLinksParam->subjRange = new CRange<TSeqPos>(from,to);
			alnLinksParam->flip = avRef->StrandSign(0) != avRef->StrandSign(1);			
		}
		else{
            TSeqPos currFrom = alnLinksParam->subjRange->GetFrom();
            TSeqPos currTo = alnLinksParam->subjRange->GetTo();
			alnLinksParam->subjRange->SetFrom(min(from,currFrom));
			alnLinksParam->subjRange->SetTo(max(to,currTo));	
		}				

	}
	if (m_AlignOption & eHtml && m_AlignOption & eShowBlastInfo) {
		alnLinksParam->hspNumber = (!first) ? alnLinksParam->hspNumber + 1 : 1;
	}
    
    if(first){		
		m_AlnLinksParams.insert(map<string, struct SAlnLinksParams>::value_type(idString,*alnLinksParam));           
	}
}

END_SCOPE(align_format)
END_NCBI_SCOPE
