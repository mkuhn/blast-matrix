/*  $Id: win_mask_config.cpp 183994 2010-02-23 20:20:11Z morgulis $
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
 *   CWinMaskConfig class member and method definitions.
 *
 */

#include <ncbi_pch.hpp>
#include <corelib/ncbidbg.hpp>

#include "win_mask_config.hpp"
#include <objtools/seqmasks_io/mask_cmdline_args.hpp>
#include <objtools/seqmasks_io/mask_fasta_reader.hpp>
#include <objtools/seqmasks_io/mask_bdb_reader.hpp>
#include <objtools/seqmasks_io/mask_writer_int.hpp>
#include <objtools/seqmasks_io/mask_writer_fasta.hpp>
#include <objtools/seqmasks_io/mask_writer_seqloc.hpp>
#include <objtools/seqmasks_io/mask_writer_blastdb_maskinfo.hpp>
#include <objects/seqloc/Seq_id.hpp>
#include <objmgr/util/sequence.hpp>

BEGIN_NCBI_SCOPE
USING_SCOPE(objects);

CMaskWriter*
CWinMaskConfig::x_GetWriter(const CArgs& args, 
                            CNcbiOstream& output, 
                            const string& format)
{
    CMaskWriter* retval = NULL;
    if (format == "interval") {
        retval = new CMaskWriterInt(output);
    } else if (format == "fasta") {
        retval = new CMaskWriterFasta(output);
    } else if (NStr::StartsWith(format, "seqloc_")) {
        retval = new CMaskWriterSeqLoc(output, format);
    } else if (NStr::StartsWith(format, "maskinfo_")) {
        retval = 
            new CMaskWriterBlastDbMaskInfo(output, format, 3,
                               eBlast_filter_program_windowmasker,
                               BuildAlgorithmParametersString(args));
    } else {
        throw runtime_error("Unknown output format");
    }
    return retval;
}

//----------------------------------------------------------------------------
CWinMaskConfig::CWinMaskConfig( const CArgs & args )
    // : is( !args["mk_counts"].AsBoolean() && args[kInputFormat].AsString() != "blastdb" ? 
    : is( !args["mk_counts"] && args[kInputFormat].AsString() != "blastdb" ? 
          // ( !args[kInput].AsString().empty() 
          ( !(args[kInput].AsString() == "-")
            ? new CNcbiIfstream( args[kInput].AsString().c_str() ) 
            : static_cast<CNcbiIstream*>(&NcbiCin) ) : NULL ), reader( NULL ), 
      // os( !args["mk_counts"].AsBoolean() ?
      os( !args["mk_counts"] ?
          // ( !args[kOutput].AsString().empty() 
          ( !(args[kOutput].AsString() == "-")
            ? new CNcbiOfstream( args[kOutput].AsString().c_str() )
            : static_cast<CNcbiOstream*>(&NcbiCout) ) : NULL ), writer( NULL ),
      // lstat_name( args["ustat"].AsString() ),
      lstat_name( args["ustat"] ? args["ustat"].AsString() : "" ),
      textend( args["t_extend"] ? args["t_extend"].AsInteger() : 0 ), 
      cutoff_score( args["t_thres"] ? args["t_thres"].AsInteger() : 0 ),
      max_score( args["t_high"] ? args["t_high"].AsInteger() : 0 ),
      min_score( args["t_low"] ? args["t_low"].AsInteger() : 0 ),
      window_size( args["window"] ? args["window"].AsInteger() : 0 ),
      // merge_pass( args["mpass"].AsBoolean() ),
      // merge_cutoff_score( args["mscore"].AsInteger() ),
      // abs_merge_cutoff_dist( args["mabs"].AsInteger() ),
      // mean_merge_cutoff_dist( args["mmean"].AsInteger() ),
      merge_pass( false ),
      merge_cutoff_score( 50 ),
      abs_merge_cutoff_dist( 8 ),
      mean_merge_cutoff_dist( 50 ),
      // trigger( args["trigger"].AsString() ),
      trigger( "mean" ),
      // tmin_count( args["tmin_count"].AsInteger() ),
      tmin_count( 0 ),
      // discontig( args["discontig"].AsBoolean() ),
      // pattern( args["pattern"].AsInteger() ),
      discontig( false ),
      pattern( 0 ),
      // window_step( args["wstep"].AsInteger() ),
      // unit_step( args["ustep"].AsInteger() ),
      window_step( 1 ),
      unit_step( 1 ),
      // merge_unit_step( args["mustep"].AsInteger() ),
      merge_unit_step( 1 ),
      // mk_counts( args["mk_counts"].AsBoolean() ),
      mk_counts( args["mk_counts"] ),
      // fa_list( args["fa_list"].AsBoolean() ),
      fa_list( args["fa_list"] ? args["fa_list"].AsBoolean() : false ),
      // mem( args["mem"].AsInteger() ),
      mem( args["mem"] ? args["mem"].AsInteger() : 1536 ),
      unit_size( args["unit"] ? args["unit"].AsInteger() : 0 ),
      genome_size( args["genome_size"] ? args["genome_size"].AsInt8() : 0 ),
      input( args[kInput].AsString() ),
      output( args[kOutput].AsString() ),
      // th( args["th"].AsString() ),
      th( "90,99,99.5,99.8" ),
      // use_dust( args["dust"].AsBoolean() ),
      use_dust( args["dust"] ? args["dust"].AsBoolean() : false ),
      // dust_window( args["dust_window"].AsInteger() ),
      // dust_linker( args["dust_linker"].AsInteger() ),
      dust_window( 64 ),
      // dust_level( 20 ),
      // dust_level( args["dust_level"].AsInteger() ),
      dust_level( args["dust_level"] ? args["dust_level"].AsInteger() : 20 ),
      dust_linker( 1 ),
      // checkdup( args["checkdup"].AsBoolean() ),
      checkdup( args["checkdup"] ? args["checkdup"].AsBoolean() : false ),
      // sformat( args["sformat"].AsString() ),
      sformat( args["sformat"] ? args["sformat"].AsString() : "ascii" ),
      // smem( args["smem"].AsInteger() ),
      smem( args["smem"] ? args["smem"].AsInteger() : 512 ),
      ids( 0 ), exclude_ids( 0 ),
      use_ba( args["use_ba"].AsBoolean() ),
      text_match( args["text_match"].AsBoolean() )
{
    _TRACE( "Entering CWinMaskConfig::CWinMaskConfig()" );

    // string iformatstr = args[kInputFormat].AsString();
    iformatstr = args[kInputFormat] ? args[kInputFormat].AsString() : "fasta";

    if( !mk_counts )
    {
        if( is && !*is )
        {
            NCBI_THROW( CWinMaskConfigException,
                        eInputOpenFail,
                        args[kInput].AsString() );
        }

        if( iformatstr == "fasta" )
            reader = new CMaskFastaReader( *is, true, args["parse_seqids"] );
        else if( iformatstr == "blastdb" )
            reader = new CMaskBDBReader( args[kInput].AsString() );

        string oformatstr = args[kOutputFormat].AsString();

        writer = x_GetWriter(args, *os, oformatstr);

        if( !reader )
        {
            NCBI_THROW( CWinMaskConfigException,
                        eReaderAllocFail, "" );
        }

        set_max_score = args["set_t_high"]  ? args["set_t_high"].AsInteger()
                                            : 0;
        set_min_score = args["set_t_low"]   ? args["set_t_low"].AsInteger()
                                            : 0;
    }
    else {
        text_match = true;
    }

    string ids_file_name( args["ids"].AsString() );
    string exclude_ids_file_name( args["exclude_ids"].AsString() );

    if(    !ids_file_name.empty()
        && !exclude_ids_file_name.empty() )
    {
        NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
                    "only one of -ids or -exclude_ids can be specified" );
    }

    if( !ids_file_name.empty() ) {
        if( text_match ) {
            ids = new CIdSet_TextMatch;
        }else {
            if( !mk_counts && iformatstr == "blastdb" ) 
                ids = new CIdSet_SeqId;
            else
                NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
                        "-text_match false can be used only with -mk_counts true "
                        "and " + string(kInputFormat) + " blastdb" );
        }

        FillIdList( ids_file_name, *ids );
    }

    if( !exclude_ids_file_name.empty() ) {
        if( text_match ) {
            exclude_ids = new CIdSet_TextMatch;
        }else {
            if( !mk_counts && iformatstr == "blastdb" ) 
                exclude_ids = new CIdSet_SeqId;
            else
                NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
                        "-text_match false can be used only with -mk_counts true "
                        "and " + string(kInputFormat) + " blastdb" );
        }

        FillIdList( exclude_ids_file_name, *exclude_ids );
    }

    _TRACE( "Leaving CWinMaskConfig::CWinMaskConfig" );
}

CWinMaskConfig::~CWinMaskConfig()
{
    if ( writer ) {
        delete writer;
    }
}

//----------------------------------------------------------------------------
void CWinMaskConfig::Validate() const
{
    _TRACE( "Entering CWinMaskConfig::Validate" );

    if( !mk_counts && lstat_name.empty() ){
        NCBI_THROW( CWinMaskConfigException, eInconsistentOptions,
                    "one of '-mk_counts true' or '-ustat <stat_file>' "
                    "must be specified" );
    }

    _TRACE( "Leaving CWinMaskConfig::Validate" );
}

//----------------------------------------------------------------------------
void CWinMaskConfig::FillIdList( const string & file_name, 
                                 CIdSet & id_list )
{
    CNcbiIfstream file( file_name.c_str() );
    string line;

    while( NcbiGetlineEOL( file, line ) ) {
        if( !line.empty() )
        {
            string::size_type stop( line.find_first_of( " \t" ) );
            string::size_type start( line[0] == '>' ? 1 : 0 );
            string id_str = line.substr( start, stop - start );
            id_list.insert( id_str );
        }
    }
}

//----------------------------------------------------------------------------
const char * CWinMaskConfig::CWinMaskConfigException::GetErrCodeString() const
{
    switch( GetErrCode() )
    {
    case eInputOpenFail: 

        return "can not open input stream";

    case eReaderAllocFail:

        return "can not allocate fasta sequence reader";

    case eInconsistentOptions:

        return "inconsistent program options";

    default: 

        return CException::GetErrCodeString();
    }
}

END_NCBI_SCOPE
