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


void PrintMatrix(Int4 **matrix, ostream& out)
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
        
        /*** Process the input ***/
        for (; !input.End(); ) {

            // Get queries
            CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
            CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));

            CRef<ILocalQueryData> local_query_data = queries->MakeLocalQueryData(&opt);
            
            BLAST_SequenceBlk* query = local_query_data->GetSequenceBlk();

            // Loop through all combinations of subjects and queries
            for (size_t query_index=0; query_index < local_query_data->GetNumQueries(); query_index++)
            {
                CConstRef<CSeq_loc> qseqloc = local_query_data->GetSeq_loc(query_index);
                
                // keep track of position within the subject BLAST_SequenceBlk
                size_t current_start = 0;
                size_t current_end = 0;
                
                // Determine query composition
                Blast_AminoAcidComposition query_composition;
                Blast_ReadAaComposition(&query_composition, BLASTAA_SIZE, query->sequence, query->length);
                
                for (size_t subject_index=0; subject_index < local_subject_data->GetNumQueries(); subject_index++)
                {
                     CConstRef<CSeq_loc> sseqloc = local_subject_data->GetSeq_loc(subject_index);
                     int subject_length = sseqloc->GetStop(eExtreme_Positional);
                     
                     current_start = current_end;
                     current_end = current_start + subject_length;
                     
                     printf("Start: %zu\n", current_start);
                     printf("End: %zu\n", current_end);

                     // Determine subject composition
                     Blast_AminoAcidComposition subject_composition;
                     Blast_ReadAaComposition(&subject_composition, BLASTAA_SIZE, &subject->sequence[current_start], subject_length);
                     
                     Blast_MatrixInfo *scaledMatrixInfo = Blast_MatrixInfoNew(BLASTAA_SIZE, BLASTAA_SIZE, 0);
                     scaledMatrixInfo->matrixName = strdup(opt.GetMatrixName());
                     scaledMatrixInfo->ungappedLambda = 0.31; // FIXME
                    
                     /* Frequency ratios for the matrix */
                     SFreqRatios * stdFreqRatios = NULL;

                     stdFreqRatios = _PSIMatrixFrequencyRatiosNew(scaledMatrixInfo->matrixName);
                     for (int i = 0;  i < BLASTAA_SIZE;  i++) {
                         for (int j = 0;  j < BLASTAA_SIZE;  j++) {
                             scaledMatrixInfo->startFreqRatios[i][j] = stdFreqRatios->data[i][j];
                         }
                     }
                     stdFreqRatios = _PSIMatrixFrequencyRatiosFree(stdFreqRatios);

                     Blast_Int4MatrixFromFreq(scaledMatrixInfo->startMatrix, scaledMatrixInfo->cols,
                                              scaledMatrixInfo->startFreqRatios, scaledMatrixInfo->ungappedLambda);
                    

                     /* which mode of composition adjustment is actually used? */
                     EMatrixAdjustRule matrix_adjust_rule = eDontAdjustMatrix;
                     double pvalueForThisPair = (-1); /* p-value for this match for composition; -1 == no adjustment*/
                     double LambdaRatio; /*lambda ratio*/

                     BlastScoreBlk *sbp = BlastScoreBlkNew(BLASTAA_SEQ_CODE, 1);
                     sbp->name = strdup(opt.GetMatrixName());
                     Blast_ScoreBlkMatrixFill(sbp, &BlastFindMatrixPath);
                     
                     PrintMatrix(sbp->matrix->data, cerr);
                     
                     int adjust_search_failed =
                         Blast_AdjustScores(sbp->matrix->data,
                                            &query_composition, query->length,
                                            &subject_composition, subject_length,
                                            scaledMatrixInfo, eCompositionBasedStats,
                                            kReMatrixAdjustmentPseudocounts, NRrecord,
                                            &matrix_adjust_rule, &s_CalcLambda,
                                            &pvalueForThisPair,
                                            compositionTestIndex,
                                            &LambdaRatio);

                    PrintMatrix(sbp->matrix->data, cerr);

                    _ASSERT(adjust_search_failed == 0);
                    
                     current_end++;
                }
            }
            

            Blast_CompositionWorkspaceFree(&NRrecord);


        }

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
