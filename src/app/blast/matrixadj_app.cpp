/*  $Id: blastp_app.cpp 194738 2010-06-16 20:14:23Z camacho $
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
 * Authors:  Christiam Camacho
 *
 */

/** @file blastp_app.cpp
 * BLASTP command line application
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
	"$Id: blastp_app.cpp 194738 2010-06-16 20:14:23Z camacho $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/remote_blast.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/blastinput/blastp_args.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/format/blast_format.hpp>
#include "blast_app_util.hpp"

#include <algo/blast/composition_adjustment/composition_adjustment.h>
#include <algo/blast/composition_adjustment/composition_constants.h>
#include "../../algo/blast/core/matrix_freq_ratios.h"
#include "../../algo/blast/api/blast_setup.hpp"

#include <objmgr/util/sequence.hpp>

#ifndef SKIP_DOXYGEN_PROCESSING
USING_NCBI_SCOPE;
USING_SCOPE(blast);
USING_SCOPE(objects);
#endif

class CBlastpApp : public CNcbiApplication
{
public:
    /** @inheritDoc */
    CBlastpApp() {
        CRef<CVersion> version(new CVersion());
        version->SetVersionInfo(new CBlastVersion());
        SetFullVersion(version);
    }
private:
    /** @inheritDoc */
    virtual void Init();
    /** @inheritDoc */
    virtual int Run();

    /// This application's command line args
    CRef<CBlastpAppArgs> m_CmdLineArgs;
};

void CBlastpApp::Init()
{
    // formulate command line arguments

    m_CmdLineArgs.Reset(new CBlastpAppArgs());

    // read the command line

    HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
    SetupArgDescriptions(m_CmdLineArgs->SetCommandLine());
}

/**
 * A callback routine: compute lambda for the given score
 * probabilities.
 * (@sa calc_lambda_type).
 */
static double
s_CalcLambda(double probs[], int min_score, int max_score, double lambda0)
{
   
    int i;                 /* loop index */      
    int score_range;       /* range of possible scores */
    double avg;            /* expected score of aligning two characters */
    Blast_ScoreFreq freq;  /* score frequency data */

    score_range = max_score - min_score + 1;
    avg = 0.0;
    for (i = 0;  i < score_range;  i++) {
        avg += (min_score + i) * probs[i];
    }
    freq.score_min = min_score;
    freq.score_max = max_score;
    freq.obs_min = min_score;
    freq.obs_max = max_score;
    freq.sprob0 = probs;
    freq.sprob = &probs[-min_score];
    freq.score_avg = avg;

    return Blast_KarlinLambdaNR(&freq, lambda0);
}

/**
 * Print a matrix on one line, replacing -32768 with X. If matrix is null, print header line
 */
void PrintMatrix(Int4 **matrix, ostream& out)
{
    for (int i = 0; i < BLASTAA_SIZE; i++)
    {
        for (int j = 0; j < BLASTAA_SIZE; j++)
        {
            if (matrix)
            {
                int m = matrix[i][j];
                if (m == -32768) {
                    out  << "\tX";
                } else {
                    out  << "\t" << m;
                }
            }
            else
            {
                out  << "\t" << NCBISTDAA_TO_AMINOACID[i] << NCBISTDAA_TO_AMINOACID[j];
            }
        }
    }
    out << "\n";
}

void PrettyPrintMatrix(Int4 **matrix, ostream& out)
{
    for (int i = 0; i < BLASTAA_SIZE; i++) out << "\t" << NCBISTDAA_TO_AMINOACID[i];
    out << "\n";

    for (int i = 0; i < BLASTAA_SIZE; i++)
    {
        out << NCBISTDAA_TO_AMINOACID[i];
        for (int j = 0; j < BLASTAA_SIZE; j++)
        {
            out  << "\t" << matrix[i][j];
        }
        out << "\n";
    }
}


void PrintSeq(Uint1* sequence, int sequence_length, ostream& out)
{
    for (int i = 0; i < sequence_length; i++) 
    {
        out << NCBISTDAA_TO_AMINOACID[ sequence[i] ];
        if (i % 80 == 79) out << "\n";
    }
    out << "\n";
}

int getCCType(const string t)
{
    const string prefix = t.substr(0, 3);
    if (prefix == "no-") return -1;
    if (prefix == "cc-") return 1;
    return 0;
}


/** Default instructions and mask residue for SEG filtering, copied from blast_kappa.c */
#define BLASTP_MASK_INSTRUCTIONS "S 10 1.8 2.1"

/**
 * Filter low complexity regions from the sequence data; uses the SEG
 * algorithm. Based on s_DoSegSequenceData in blast_kappa.c
 *
 * @param seqData            data to be filtered
 * @return   0 for success; -1 for out-of-memory
 */
static int
MaskSequence(Uint1* sequence, int sequence_length)
{
    int status = 0;
    BlastSeqLoc* mask_seqloc = NULL;
    SBlastFilterOptions* filter_options = NULL;

    status = BlastFilteringOptionsFromString(eBlastTypeBlastp,
                                             BLASTP_MASK_INSTRUCTIONS,
                                             &filter_options, NULL);
    if (status == 0) {
        status = BlastSetUp_Filter(eBlastTypeBlastp, sequence,
                                   sequence_length, 0, filter_options,
                                   &mask_seqloc, NULL);
        filter_options = SBlastFilterOptionsFree(filter_options);
    }
    if (status == 0) {
        Blast_MaskTheResidues(sequence, sequence_length,
                              FALSE, mask_seqloc, FALSE, 0);
    }
    if (mask_seqloc != NULL) {
        mask_seqloc = BlastSeqLocFree(mask_seqloc);
    }
    return status;
}


int CBlastpApp::Run(void)
{
    int status = BLAST_EXIT_SUCCESS;

    try {

        // Allow the fasta reader to complain on invalid sequence input
        SetDiagPostLevel(eDiag_Warning);

        /*** Get the BLAST options ***/
        const CArgs& args = GetArgs();
        RecoverSearchStrategy(args, m_CmdLineArgs);
        CRef<CBlastOptionsHandle> opts_hndl(&*m_CmdLineArgs->SetOptions(args));
        const CBlastOptions& opt = opts_hndl->GetOptions();

        /*** Get the query sequence(s) ***/
        CRef<CQueryOptionsArgs> query_opts = 
            m_CmdLineArgs->GetQueryOptionsArgs();
        SDataLoaderConfig dlconfig(query_opts->QueryIsProtein());
        dlconfig.OptimizeForWholeLargeSequenceRetrieval();
        CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
                                     query_opts->UseLowercaseMasks(),
                                     query_opts->GetParseDeflines(),
                                     query_opts->GetRange(),
                                     !m_CmdLineArgs->ExecuteRemotely());
        iconfig.SetQueryLocalIdMode();
        CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
        CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());

        /*** Initialize the database/subject ***/
        CRef<CBlastDatabaseArgs> db_args(m_CmdLineArgs->GetBlastDatabaseArgs());
        CRef<CLocalDbAdapter> db_adapter;
        CRef<CScope> scope;
        InitializeSubject(db_args, opts_hndl, m_CmdLineArgs->ExecuteRemotely(),
                         db_adapter, scope);
        _ASSERT(db_adapter && scope);

        // Get a IQueryFactory from our CLocalDbAdapter, using a custom method added to this class
        CRef<IQueryFactory> subjects = db_adapter->GetSubjectFactory();
        CRef<ILocalQueryData> local_subject_data = subjects->MakeLocalQueryData(&opt);
        BLAST_SequenceBlk* subject = local_subject_data->GetSequenceBlk();
        
        // Setup for composition-based adjustment
        Blast_CompositionWorkspace *NRrecord = Blast_CompositionWorkspaceNew();
        int status_code = Blast_CompositionWorkspaceInit(NRrecord, opt.GetMatrixName());

        _ASSERT(status_code == 0);
        
        static const int kReMatrixAdjustmentPseudocounts = 20;
        static const int compositionTestIndex = 0;
        
        static const int scaling_factor = 32;
        
        // print header
        fflush(stdout);
        cout << "# query\tsubject\tscaling_factor";
        PrintMatrix(0, cout);
        cout.flush();
        
        /*** Process the input ***/
        for (; !input.End(); ) {

            // Get queries
            CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
            CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));

            CRef<ILocalQueryData> local_query_data = queries->MakeLocalQueryData(&opt);
            
            BLAST_SequenceBlk* query = local_query_data->GetSequenceBlk();

            // keep track of position within the query BLAST_SequenceBlk
            size_t query_start = 0;
            size_t query_end = 0;

            // Loop through all combinations of subjects and queries
            for (size_t query_index=0; query_index < local_query_data->GetNumQueries(); query_index++)
            {
                CConstRef<CSeq_loc> qseqloc = local_query_data->GetSeq_loc(query_index);
                
                const string query_title = sequence::GetTitle( scope->GetBioseqHandle(*qseqloc) );
                const int query_type = getCCType(query_title);
                
                const int query_length = 1 + qseqloc->GetStop(eExtreme_Positional);
                query_start = query_end;
                query_end = query_start + query_length;
                
                // Determine query composition
                Blast_AminoAcidComposition query_composition;
                Blast_ReadAaComposition(&query_composition, BLASTAA_SIZE, &query->sequence[query_start], query_length);
                
                // fflush(stdout);
                // cout << "query:\n";
                // PrintSeq(&query->sequence[query_start], query_length, cout);
                // cout << "query_start: " << query_start << " query_length " << query_length << "\n";
                // cout.flush();
                
                // keep track of position within the subject BLAST_SequenceBlk
                size_t subject_start = 0;
                size_t subject_end = 0;
                
                for (size_t subject_index=0; subject_index < local_subject_data->GetNumQueries(); subject_index++)
                {
                     CConstRef<CSeq_loc> sseqloc = local_subject_data->GetSeq_loc(subject_index);
                     
                     string subject_title = sequence::GetTitle( scope->GetBioseqHandle(*sseqloc) );

                     const int subject_length = 1 + sseqloc->GetStop(eExtreme_Positional);
                     subject_start = subject_end;
                     subject_end = subject_start + subject_length;
                     
                     // if a query/subject type is set (i.e. prefix "cc-" / "no-"), only calculate matrices for same type
                     const int subject_type = getCCType(subject_title);
                     if (subject_type == query_type) 
                     {
                          MaskSequence(&subject->sequence[subject_start], subject_length);

                          BlastScoreBlk *sbp = BlastScoreBlkNew(BLASTAA_SEQ_CODE, 1);
                          sbp->name = strdup(opt.GetMatrixName());
                          Blast_ScoreBlkMatrixFill(sbp, &BlastFindMatrixPath);
                          Blast_ScoreBlkKbpIdealCalc(sbp);

                          // Determine subject composition
                          Blast_AminoAcidComposition subject_composition;
                          Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, &subject->sequence[subject_start], subject_length);

                          // fflush(stdout);
                          // cout << "subject:\n";
                          // PrintSeq(&subject->sequence[subject_start], subject_length, cout);
                          // cout << "subject_start: " << subject_start << " subject_length " << subject_length << "\n";
                          // cout.flush();

                          Blast_MatrixInfo *scaledMatrixInfo = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, 0);
                          scaledMatrixInfo->matrixName = strdup(opt.GetMatrixName());
                          scaledMatrixInfo->ungappedLambda = sbp->kbp_ideal->Lambda / scaling_factor;

                          /* Frequency ratios for the matrix */
                          SFreqRatios * stdFreqRatios = NULL;

                          stdFreqRatios = _PSIMatrixFrequencyRatiosNew(scaledMatrixInfo->matrixName);
                          for (int i = 0;  i < BLASTAA_SIZE;  i++) {
                              for (int j = 0;  j < BLASTAA_SIZE;  j++) {
                                  scaledMatrixInfo->startFreqRatios[i][j] = stdFreqRatios->data[i][j];
                              }
                          }
                          _PSIMatrixFrequencyRatiosFree(stdFreqRatios);

                          Blast_Int4MatrixFromFreq(scaledMatrixInfo->startMatrix, scaledMatrixInfo->cols,
                                                   scaledMatrixInfo->startFreqRatios, scaledMatrixInfo->ungappedLambda);


                          /* which mode of composition adjustment is actually used? */
                          EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;
                          double pvalueForThisPair = (-1); /* p-value for this match for composition; -1 == no adjustment*/
                          double LambdaRatio; /*lambda ratio*/

                          int adjust_search_failed =
                              Blast_AdjustScores(sbp->matrix->data,
                                                 &query_composition, query_length,
                                                 &subject_composition, subject_length,
                                                 scaledMatrixInfo, eCompositionMatrixAdjust,
                                                 kReMatrixAdjustmentPseudocounts, NRrecord,
                                                 &matrix_adjust_rule, &s_CalcLambda,
                                                 &pvalueForThisPair,
                                                 compositionTestIndex,
                                                 &LambdaRatio);

                         // PrettyPrintMatrix(sbp->matrix->data, cout);

                         // if the adjustment failed, don't print anything: we'll fall back to other matrices in CCAlign
                         if (!adjust_search_failed)
                         {
                             fflush(stdout);
                             cout << query_title << "\t" << subject_title << "\t" << scaling_factor;
                             PrintMatrix(sbp->matrix->data, cout);
                             cout.flush();                         
                         }
                         
                         Blast_MatrixInfoFree(&scaledMatrixInfo);
                         BlastScoreBlkFree(sbp);
                     }
                     
                     subject_end++;
                }
                query_end++;
            }
        }

        Blast_CompositionWorkspaceFree(&NRrecord);

        if (m_CmdLineArgs->ProduceDebugOutput()) {
            opts_hndl->GetOptions().DebugDumpText(NcbiCerr, "BLAST options", 1);
        }

    } CATCH_ALL(status)
    return status;
}

#ifndef SKIP_DOXYGEN_PROCESSING
int main(int argc, const char* argv[] /*, const char* envp[]*/)
{
    return CBlastpApp().AppMain(argc, argv, 0, eDS_Default, 0);
}
#endif /* SKIP_DOXYGEN_PROCESSING */
