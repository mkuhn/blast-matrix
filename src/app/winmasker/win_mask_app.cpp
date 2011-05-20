/*  $Id: win_mask_app.cpp 183994 2010-02-23 20:20:11Z morgulis $
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
 * Author:  Aleksandr Morgulis
 *
 * File Description:
 *   CWinMaskApplication class member and method definitions.
 *
 */

#include <ncbi_pch.hpp>
#include <corelib/ncbidbg.hpp>
#include <objtools/readers/fasta.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Seq_inst.hpp>
#include <objects/seq/Seq_data.hpp>
#include <objects/seq/seqport_util.hpp>
#include <objects/seq/IUPACna.hpp>
#include <objects/seqloc/Seq_id.hpp>

#include <objmgr/bioseq_ci.hpp>
#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/seq_entry_handle.hpp>
#include <objtools/data_loaders/genbank/gbloader.hpp>

#include <objtools/seqmasks_io/mask_cmdline_args.hpp>
#include <objtools/seqmasks_io/mask_reader.hpp>
#include <objtools/seqmasks_io/mask_fasta_reader.hpp>
#include <objtools/seqmasks_io/mask_writer.hpp>
#include <algo/winmask/seq_masker.hpp>
#include <algo/winmask/win_mask_gen_counts.hpp>
#include <algo/winmask/win_mask_util.hpp>
#include <algo/winmask/win_mask_counts_converter.hpp>
#include "win_mask_app.hpp"
#include "win_mask_config.hpp"
#include "win_mask_sdust_masker.hpp"


BEGIN_NCBI_SCOPE
USING_SCOPE(objects);

//-------------------------------------------------------------------------
const char * const 
CWinMaskApplication::USAGE_LINE = "Window based sequence masker";

//-------------------------------------------------------------------------
void CWinMaskApplication::Init(void)
{
    HideStdArgs(fHideLogfile | fHideConffile | fHideVersion | fHideDryRun);
    auto_ptr< CArgDescriptions > arg_desc( new CArgDescriptions );

    // Set the program description
    arg_desc->SetUsageContext( GetArguments().GetProgramBasename(),
                               USAGE_LINE );

    // Adding command line arguments descriptions
    arg_desc->AddDefaultKey( "ustat", "unit_counts",
                             "file with unit counts",
                             CArgDescriptions::eString, "" );
    /*
    arg_desc->AddDefaultKey( kInput, "input_file_name",
                             "input file name "
                             "(not optional if used with -mk_counts option)",
                             CArgDescriptions::eString, "" );
    arg_desc->AddDefaultKey( kOutput, "output_file_name",
                             "output file name",
                             CArgDescriptions::eString, "" );
    */
    arg_desc->AddDefaultKey( kInput, "input_file_name",
                             "input file name "
                             "(not optional if used with -mk_counts or -convert options)",
                             CArgDescriptions::eInputFile, "-" );
    arg_desc->AddDefaultKey( kOutput, "output_file_name",
                             "output file name",
                             CArgDescriptions::eOutputFile, "-" );
    arg_desc->AddDefaultKey( "checkdup", "check_duplicates",
                             "check for duplicate sequences",
                             CArgDescriptions::eBoolean, "false" );
    arg_desc->AddOptionalKey( "window", "window_size", "window size",
                              CArgDescriptions::eInteger );
#if 0
    arg_desc->AddDefaultKey( "wstep", "window_step", "window step",
                             CArgDescriptions::eInteger, "1" );
    arg_desc->AddDefaultKey( "ustep", "unit_step", "unit step",
                             CArgDescriptions::eInteger, "1" );
#endif
    arg_desc->AddOptionalKey( "t_extend", "T_extend", 
                              "window score above which it is allowed to extend masking",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "t_thres", "T_threshold",
                              "window score threshold used to trigger masking",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "t_high", "T_high",
                              "maximum useful unit score",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "t_low", "T_low",
                              "minimum useful unit score",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "set_t_high", "score_value",
                              "alternative high score for a unit if the"
                              "original unit score is more than highscore",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "set_t_low", "score_value",
                              "alternative low score for a unit if the"
                              "original unit score is lower than lowscore",
                              CArgDescriptions::eInteger );
#if 0
    arg_desc->AddDefaultKey( "ambig", "ambiguity_handler",
                             "the way to handle ambiguity characters",
                             CArgDescriptions::eString, "break" );
#endif
    arg_desc->AddDefaultKey( kInputFormat, "input_format",
                             "controls the format of the masker input (for masking stage only)",
                             CArgDescriptions::eString, *kInputFormats );
    arg_desc->AddFlag      ( "parse_seqids",
                             "Parse Seq-ids in FASTA input", true );
    arg_desc->AddDefaultKey( kOutputFormat, "output_format",
                             "controls the format of the masker output (for masking stage only)",
                             CArgDescriptions::eString, *kOutputFormats );
    arg_desc->AddDefaultKey( "sformat", "unit_counts_format",
                             "controls the format of the output file containing the unit counts "
                             "(for counts generation and conversion only)",
                             CArgDescriptions::eString, "ascii" );
#if 0
    arg_desc->AddDefaultKey( "mpass", "merge_pass_flag",
                             "true if separate merging pass is needed",
                             CArgDescriptions::eBoolean, "false" );
    arg_desc->AddDefaultKey( "discontig", "discontiguous_units",
                             "true if using discontiguous units",
                             CArgDescriptions::eBoolean, "false" );
    arg_desc->AddDefaultKey( "mscore", "merge_cutoff_score",
                             "minimum average unit score triggering a merge",
                             CArgDescriptions::eInteger, "50" );
    arg_desc->AddDefaultKey( "mabs", "distance",
                             "absolute distance threshold for merging",
                             CArgDescriptions::eInteger, "8" );
    arg_desc->AddDefaultKey( "mmean", "distance",
                             "distance threshold for merging if average unit"
                             " score is high enough",
                             CArgDescriptions::eInteger, "50" );
    arg_desc->AddDefaultKey( "mustep", "merge_unit_step",
                             "unit step value used for interval merging",
                             CArgDescriptions::eInteger, "1" );
    arg_desc->AddDefaultKey( "trigger", "trigger_type",
                             "type of the event triggering masking",
                             CArgDescriptions::eString, "mean" );
    arg_desc->AddDefaultKey( "tmin_count", "unit_count",
                             "number of units to count with min trigger",
                             CArgDescriptions::eInteger, "0" );
    arg_desc->AddDefaultKey( "pattern", "base_mask", 
                             "which bases in a window to use as a discontinuous unit",
                             CArgDescriptions::eInteger, "0" );
    arg_desc->AddDefaultKey( "dbg", "debug_output",
                             "enable debug output",
                             CArgDescriptions::eBoolean, "false" );
#endif
    /*
    arg_desc->AddDefaultKey( "mk_counts", "generate_counts",
                             "generate frequency counts for a database",
                             CArgDescriptions::eBoolean, "false" );
    */
    arg_desc->AddFlag( "mk_counts", "generate frequency counts for a database" );
    arg_desc->AddFlag( "convert", "convert counts between different formats" );

    /*
    arg_desc->AddDefaultKey( "convert", "convert_counts",
                             "convert counts between different formats",
                             CArgDescriptions::eBoolean, "false" );
    */
    arg_desc->AddDefaultKey( "fa_list", "input_is_a_list",
                             "indicates that -input represents a file containing "
                             "a list of names of fasta files to process, one name "
                             " per line", 
                             CArgDescriptions::eBoolean, "false" );
    arg_desc->AddDefaultKey( "mem", "available_memory",
                             "memory available for mk_counts option in megabytes",
                             CArgDescriptions::eInteger, "1536" );
    arg_desc->AddDefaultKey( "smem", "available_memory",
                             "memory available for masking stage in megabytes",
                             CArgDescriptions::eInteger, "512" );
    arg_desc->AddOptionalKey( "unit", "unit_length",
                              "number of bases in a unit",
                              CArgDescriptions::eInteger );
    arg_desc->AddOptionalKey( "genome_size", "genome_size",
                              "total size of the genome",
                              CArgDescriptions::eInteger );
#if 0
    arg_desc->AddDefaultKey( "th", "thresholds",
                             "4 percentage values used to determine "
                             "masking thresholds (4 floating point numbers "
                             "separated by commas)", 
                             CArgDescriptions::eString, "90,99,99.5,99.8" );
#endif
    arg_desc->AddDefaultKey( "dust", "use_dust",
                             "combine window masking with dusting",
                             CArgDescriptions::eBoolean, "F" );
    arg_desc->AddDefaultKey( "dust_level", "dust_level",
                             "dust minimum level",
                             CArgDescriptions::eInteger, "20" );
#if 0
    arg_desc->AddDefaultKey( "dust_window", "dust_window",
                             "window size for dusting",
                             CArgDescriptions::eInteger, "64" );
    arg_desc->AddDefaultKey( "dust_linker", "dust_linker",
                             "link windows by this many basepairs",
                             CArgDescriptions::eInteger, "1" );
#endif
    arg_desc->AddDefaultKey( "exclude_ids", "exclude_id_list",
                             "file containing the list of ids to exclude from processing",
                             CArgDescriptions::eString, "" );
    arg_desc->AddDefaultKey( "ids", "id_list",
                             "file containing the list of ids to process",
                             CArgDescriptions::eString, "" );
    arg_desc->AddDefaultKey( "text_match", "text_match_ids",
                             "match ids as strings",
                             CArgDescriptions::eBoolean, "T" );
    arg_desc->AddDefaultKey( "use_ba", "use_bit_array_optimization",
                             "whether to use extra bit array optimization "
                             "for optimized binary counts format",
                             CArgDescriptions::eBoolean, "T" );

    // Set some constraints on command line parameters
    arg_desc->SetConstraint( "window",
                             new CArgAllow_Integers( 1, kMax_Int ) );
#if 0
    arg_desc->SetConstraint( "wstep",
                             new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "ustep",
                             new CArgAllow_Integers( 1, 256 ) );
#endif
    arg_desc->SetConstraint( "t_extend",
                             new CArgAllow_Integers( 0, kMax_Int ) );
    arg_desc->SetConstraint( "t_thres",
                             new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "t_high",
                             new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "t_low",
                             new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "set_t_high",
                             new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "set_t_low",
                             new CArgAllow_Integers( 1, kMax_Int ) );
#if 0
    arg_desc->SetConstraint( "mscore",
                             new CArgAllow_Integers( 0, kMax_Int ) );
    arg_desc->SetConstraint( "mabs",
                             new CArgAllow_Integers( 0, kMax_Int ) );
    arg_desc->SetConstraint( "mmean",
                             new CArgAllow_Integers( 0, kMax_Int ) );
    arg_desc->SetConstraint( "mustep",
                             new CArgAllow_Integers( 0, 256 ) );
    arg_desc->SetConstraint( "ambig",
                             (new CArgAllow_Strings())->Allow( "break" ) );
#endif
    CArgAllow_Strings* strings_allowed = new CArgAllow_Strings();
    for (size_t i = 0; i < kNumInputFormats; i++) {
        strings_allowed->Allow(kInputFormats[i]);
    }
    arg_desc->SetConstraint( kInputFormat, strings_allowed );
    strings_allowed = new CArgAllow_Strings();
    for (size_t i = 0; i < kNumOutputFormats; i++) {
        strings_allowed->Allow(kOutputFormats[i]);
    }
    arg_desc->SetConstraint( kOutputFormat, strings_allowed );
    arg_desc->SetConstraint( "sformat",
                             (new CArgAllow_Strings())
                             ->Allow( "ascii" )
                             ->Allow( "binary" )
                             ->Allow( "oascii" )
                             ->Allow( "obinary" ) );
#if 0
    arg_desc->SetConstraint( "trigger",
                             (new CArgAllow_Strings())->Allow( "mean" )
                             ->Allow( "min" ) );
    arg_desc->SetConstraint( "tmin_count",
                             new CArgAllow_Integers( 0, kMax_Int ) );
#endif
    arg_desc->SetConstraint( "mem", new CArgAllow_Integers( 1, kMax_Int ) );
    arg_desc->SetConstraint( "unit", new CArgAllow_Integers( 1, 16 ) );

    // Set up dependencies.
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "outfmt" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "ustat" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "window" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "t_thres" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "t_extend" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "set_t_low" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "set_t_high" );
    // arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "infmt" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "dust" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "dust_level" );
    arg_desc->CArgDescriptions::SetDependency( "mk_counts", CArgDescriptions::eExcludes, "convert" );

    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "checkdup" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "fa_list" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "mem" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "unit" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "genome_size" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "sformat" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "smem" );
    arg_desc->CArgDescriptions::SetDependency( "ustat", CArgDescriptions::eExcludes, "convert" );

    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "checkdup" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "window" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "t_extend" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "t_thres" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "t_high" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "t_low" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "set_t_low" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "set_t_high" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "infmt" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "outfmt" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "parse_seqids" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "fa_list" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "mem" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "unit" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "genome_size" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "dust" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "dust_level" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "exclude_ids" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "ids" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "text_match" );
    arg_desc->CArgDescriptions::SetDependency( "convert", CArgDescriptions::eExcludes, "use_ba" );

    // Parse the arguments according to descriptions.
    SetupArgDescriptions(arg_desc.release());
}

//-------------------------------------------------------------------------
int CWinMaskApplication::Run (void)
{
    SetDiagPostLevel( eDiag_Warning );

    // Branch away immediately if the converter is called.
    //
    // if( GetArgs()["convert"].AsBoolean() ) {
    if( GetArgs()["convert"] ) {
        if( GetArgs()[kOutput].AsString() == "-" ) {
            CWinMaskCountsConverter converter( 
                    GetArgs()[kInput].AsString(),
                    NcbiCout,
                    GetArgs()["sformat"].AsString() + 
                    GetArgs()["smem"].AsString() );
            return converter();
        }
        else {
            CWinMaskCountsConverter converter( 
                    GetArgs()[kInput].AsString(),
                    GetArgs()[kOutput].AsString(),
                    GetArgs()["sformat"].AsString() + 
                    GetArgs()["smem"].AsString() );
            return converter();
        }
    }

    CRef<CObjectManager> om(CObjectManager::GetInstance());

#if 0
    CGBDataLoader::RegisterInObjectManager(
        *om, "id2", CObjectManager::eDefault );

    if( GetArgs()["dbg"].AsBoolean() )
        SetDiagTrace( eDT_Enable );
#endif

    // Read and validate configuration values.
    CWinMaskConfig aConfig( GetArgs() );
    aConfig.Validate();

    if( aConfig.MakeCounts() )
    {
        if( aConfig.Output() == "-" ) {
            CWinMaskCountsGenerator cg( aConfig.Input(),
                                        NcbiCout,
                                        aConfig.InFmt(),
                                        aConfig.SFormat(),
                                        aConfig.Th(),
                                        aConfig.Mem(),
                                        aConfig.UnitSize(),
                                        aConfig.GenomeSize(),
                                        aConfig.MinScore(),
                                        aConfig.MaxScore(),
                                        aConfig.CheckDup(),
                                        aConfig.FaList(),
                                        aConfig.Ids(),
                                        aConfig.ExcludeIds(),
                                        aConfig.UseBA() );
            cg();
        }
        else {
            CWinMaskCountsGenerator cg( aConfig.Input(),
                                        aConfig.Output(),
                                        aConfig.InFmt(),
                                        aConfig.SFormat(),
                                        aConfig.Th(),
                                        aConfig.Mem(),
                                        aConfig.UnitSize(),
                                        aConfig.GenomeSize(),
                                        aConfig.MinScore(),
                                        aConfig.MaxScore(),
                                        aConfig.CheckDup(),
                                        aConfig.FaList(),
                                        aConfig.Ids(),
                                        aConfig.ExcludeIds(),
                                        aConfig.UseBA() );
            cg();
        }

        return 0;
    }

    CMaskReader & theReader = aConfig.Reader();
    CMaskWriter & theWriter = aConfig.Writer();
    CSeqMasker theMasker( aConfig.LStatName(),
                          aConfig.WindowSize(),
                          aConfig.WindowStep(),
                          aConfig.UnitStep(),
                          aConfig.Textend(),
                          aConfig.CutoffScore(),
                          aConfig.MaxScore(),
                          aConfig.MinScore(),
                          aConfig.SetMaxScore(),
                          aConfig.SetMinScore(),
                          aConfig.MergePass(),
                          aConfig.MergeCutoffScore(),
                          aConfig.AbsMergeCutoffDist(),
                          aConfig.MeanMergeCutoffDist(),
                          aConfig.MergeUnitStep(),
                          aConfig.Trigger(),
                          aConfig.TMin_Count(),
                          aConfig.Discontig(),
                          aConfig.Pattern(),
                          aConfig.UseBA() );
    CRef< CSeq_entry > aSeqEntry( 0 );
    Uint4 total = 0, total_masked = 0;
    CSDustMasker * duster( 0 );
    const CWinMaskConfig::CIdSet * ids( aConfig.Ids() );
    const CWinMaskConfig::CIdSet * exclude_ids( aConfig.ExcludeIds() );

    if( aConfig.UseDust() )
        duster = new CSDustMasker( aConfig.DustWindow(),
                                   aConfig.DustLevel(),
                                   aConfig.DustLinker() );

    while( (aSeqEntry = theReader.GetNextSequence()).NotEmpty() )
    {
        if( aSeqEntry->Which() == CSeq_entry::e_not_set ) continue;
        CScope scope(*om);
        CSeq_entry_Handle seh = scope.AddTopLevelSeqEntry(*aSeqEntry);
        Uint4 masked = 0;
        CBioseq_CI bs_iter(seh, CSeq_inst::eMol_na);
        for ( ;  bs_iter;  ++bs_iter) {
            CBioseq_Handle bsh = *bs_iter;
            if (bsh.GetBioseqLength() == 0) {
                continue;
            }

            if( CWinMaskUtil::consider( bsh, ids, exclude_ids ) )
            {
                TSeqPos len = bsh.GetBioseqLength();
                total += len;
                _TRACE( "Sequence length " << len );
                CSeqVector data =
                    bsh.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
                auto_ptr< CSeqMasker::TMaskList > mask_info( theMasker( data ) );
                CSeqMasker::TMaskList dummy;

                if( duster != 0 ) // Dust and merge with mask_info
                {
                    auto_ptr< CSeqMasker::TMaskList > dust_info( 
                        (*duster)( data, *mask_info.get() ) );
                    CSeqMasker::MergeMaskInfo( mask_info.get(), dust_info.get() );
                }

                // theWriter.Print( bsh, *mask_info, aConfig.MatchId() );
                theWriter.Print( bsh, *mask_info, GetArgs()["parse_seqids"] );

                for( CSeqMasker::TMaskList::const_iterator i = mask_info->begin();
                     i != mask_info->end(); ++i )
                    masked += i->second - i->first + 1;

                total_masked += masked;
                _TRACE( "Number of positions masked: " << masked );
            }
        }
    }

    _TRACE( "Total number of positions: " << total );
    _TRACE( "Total number of positions masked: " << total_masked );
    return 0;
}

END_NCBI_SCOPE
