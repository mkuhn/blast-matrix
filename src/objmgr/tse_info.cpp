/*  $Id: tse_info.cpp 205152 2010-09-15 14:09:19Z vasilche $
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
*   TSE info -- entry for data source seq-id to TSE map
*
*/


#include <ncbi_pch.hpp>
#include <objmgr/impl/data_source.hpp>
#include <objmgr/impl/tse_info.hpp>
#include <objmgr/impl/tse_split_info.hpp>
#include <objmgr/impl/tse_chunk_info.hpp>
#include <objmgr/impl/tse_assigner.hpp>
#include <objmgr/impl/bioseq_info.hpp>
#include <objmgr/impl/bioseq_set_info.hpp>
#include <objmgr/impl/annot_object.hpp>
#include <objmgr/impl/seq_annot_info.hpp>
#include <objmgr/impl/snp_annot_info.hpp>
#include <objmgr/impl/annot_type_index.hpp>
#include <objmgr/impl/handle_range.hpp>

#include <objects/seqset/Seq_entry.hpp>

#include <objmgr/objmgr_exception.hpp>
#include <objmgr/error_codes.hpp>

#include <algorithm>


#define NCBI_USE_ERRCODE_X   ObjMgr_TSEinfo

BEGIN_NCBI_SCOPE

NCBI_DEFINE_ERR_SUBCODE_X(3);

BEGIN_SCOPE(objects)


SIdAnnotObjs::SIdAnnotObjs(void)
{
}


SIdAnnotObjs::~SIdAnnotObjs(void)
{
    NON_CONST_ITERATE ( TAnnotSet, it, m_AnnotSet ) {
        delete *it;
        *it = 0;
    }
}


SIdAnnotObjs::TRangeMap& SIdAnnotObjs::x_GetRangeMap(size_t index)
{
    if ( index >= m_AnnotSet.size() ) {
        m_AnnotSet.resize(index+1);
    }
    TRangeMap*& slot = m_AnnotSet[index];
    if ( !slot ) {
        slot = new TRangeMap;
    }
    return *slot;
}


bool SIdAnnotObjs::x_CleanRangeMaps(void)
{
    while ( !m_AnnotSet.empty() ) {
        TRangeMap*& slot = m_AnnotSet.back();
        if ( slot ) {
            if ( !slot->empty() ) {
                return false;
            }
            delete slot;
            slot = 0;
        }
        m_AnnotSet.pop_back();
    }
    return true;
}


SIdAnnotObjs::SIdAnnotObjs(const SIdAnnotObjs& _DEBUG_ARG(objs))
{
    _ASSERT(objs.m_AnnotSet.empty());
    _ASSERT(objs.m_SNPSet.empty());
}


////////////////////////////////////////////////////////////////////
//
//  CTSE_Info::
//
//    General information and indexes for top level seq-entries
//


CTSE_Info::CTSE_Info(void) 
    : m_InternalBioObjNumber(0)
{
    x_Initialize();

    x_TSEAttach(*this);
}


CTSE_Info::CTSE_Info(const TBlobId& blob_id,
                     TBlobVersion blob_version)
    : m_InternalBioObjNumber(0)

{
    x_Initialize();

    m_BlobId = blob_id;
    m_BlobVersion = blob_version;

    x_TSEAttach(*this);
}


CTSE_Info::CTSE_Info(CSeq_entry& entry,
                     TBlobState blob_state,
                     const TBlobId& blob_id,
                     TBlobVersion blob_version)
    : m_InternalBioObjNumber(0)

{
    x_Initialize();

    m_BlobId = blob_id;
    m_BlobVersion = blob_version;
    m_BlobState = blob_state;

    SetSeq_entry(entry);
    m_LoadState = eLoaded;

    x_TSEAttach(*this);
}


CTSE_Info::CTSE_Info(CSeq_entry& entry,
                     TBlobState blob_state)
    : m_InternalBioObjNumber(0)

{
    x_Initialize();

    m_BlobState = blob_state;

    SetSeq_entry(entry);
    m_LoadState = eLoaded;

    x_TSEAttach(*this);
}


CTSE_Info::CTSE_Info(const CTSE_Lock& tse)
    : m_BaseTSE(new SBaseTSE(tse)),
      m_InternalBioObjNumber(0)

{
    x_Initialize();

    m_BlobState = tse->m_BlobState;
    m_Name = tse->m_Name;
    m_UsedMemory = tse->m_UsedMemory;
    m_LoadState = eLoaded;


    x_SetObject(*tse, &m_BaseTSE->m_ObjectCopyMap);
    x_TSEAttach(*this);

    m_Split = tse->m_Split;
    if (m_Split) {
        CRef<ITSE_Assigner> lsnr = m_Split->GetAssigner(*tse);
        if( !lsnr )
            lsnr.Reset(new CTSE_Default_Assigner);
        m_Split->x_TSEAttach(*this, lsnr);
    }
    if (tse->HasDataSource()) {
        CDataLoader* ld = tse->GetDataSource().GetDataLoader();
        if ( ld ) {
            m_EditSaver = ld->GetEditSaver();
            m_BlobId = tse->m_BlobId;
            //CBlobIdKey(new CBlobIdString(tse->m_BlobId.ToString()));
        }
    }
}


CTSE_Info::~CTSE_Info(void)
{   
    _ASSERT(m_LockCounter.Get() == 0);
    _ASSERT(m_DataSource == 0);
    if( m_Split )
        m_Split->x_TSEDetach(*this);
}


CTSE_Info& CTSE_Info::Assign(const CTSE_Lock& tse)
{
    //    m_BaseTSE.reset(new SBaseTSE(tse));
    m_BlobState = tse->m_BlobState;
    m_Name = tse->m_Name;
    m_UsedMemory = tse->m_UsedMemory;

    if (tse->m_Contents)
        x_SetObject(*tse,NULL);//tse->m_BaseTSE->m_ObjectCopyMap);
    //x_TSEAttach(*this);

    m_Split = tse->m_Split;
    if (m_Split) {
        CRef<ITSE_Assigner> listener = m_Split->GetAssigner(*tse);
        if( !listener ) 
            listener.Reset(new CTSE_Default_Assigner);
        m_Split->x_TSEAttach(*this, listener);
    }
    return *this;
}


CTSE_Info& CTSE_Info::Assign(const CTSE_Lock& tse, 
                             CRef<CSeq_entry> entry)
//                             CRef<ITSE_Assigner> listener)
{
    m_BlobState = tse->m_BlobState;
    m_Name = tse->m_Name;
    m_UsedMemory = tse->m_UsedMemory;

    if (entry)
        SetSeq_entry(*entry);

    m_Split = tse->m_Split;
    if (m_Split) {
        CRef<ITSE_Assigner> listener = m_Split->GetAssigner(*tse);
        if( !listener ) 
            listener.Reset(new CTSE_Default_Assigner);
        m_Split->x_TSEAttach(*this, listener);
        /*
        if( !listener ) {
            listener = m_Split->GetAssigner(*tse);
            if( !listener ) {
                listener.Reset(new CTSE_Default_Assigner);
            }
        }
        m_Split->x_TSEAttach(*this, listener);
        */
    }

    //x_TSEAttach(*this);
    return *this;
}


void CTSE_Info::x_Initialize(void)
{
    m_DataSource = 0;
    m_BlobVersion = -1;
    m_BlobState = CBioseq_Handle::fState_none;
    m_UsedMemory = 0;
    m_LoadState = eNotLoaded;
    m_CacheState = eNotInCache;
    m_AnnotIdsFlags = 0;
}


void CTSE_Info::SetBlobVersion(TBlobVersion version)
{
    _ASSERT(version >= 0);
    _ASSERT(m_LoadState == eNotLoaded || !m_Object ||
            m_BlobVersion < 0 || m_BlobVersion == version);
    m_BlobVersion = version;
}


void CTSE_Info::SetName(const CAnnotName& name)
{
    m_Name = name;
}


void CTSE_Info::SetUsedMemory(size_t size)
{
    m_UsedMemory = size;
}


void CTSE_Info::SetSeq_entry(CSeq_entry& entry, CTSE_SetObjectInfo* set_info)
{
    if ( Which() != CSeq_entry::e_not_set ) {
        if ( m_LoadState == eNotLoaded ) {
            Reset();
            m_Object.Reset();
            m_Split.Reset();
            m_RequestedId.Reset();
            m_Removed_Bioseq_sets.clear();
            m_Removed_Bioseqs.clear();
            m_AnnotIdsFlags = 0;
        }
    }

    entry.Parentize();

    m_SetObjectInfo = set_info;
    if ( HasDataSource() ) {
        {{
            CDataSource::TMainLock::TWriteLockGuard guard
                (GetDataSource().m_DSMainLock);
            x_SetObject(entry);
        }}
        UpdateAnnotIndex();
    }
    else {
        x_SetObject(entry);
    }
    if ( set_info ) {
        if ( !set_info->m_Seq_annot_InfoMap.empty() ) {
            NCBI_THROW(CObjMgrException, eAddDataError,
                       "Unknown SNP annots");
        }
        m_SetObjectInfo = null;
    }
}


CBioObjectId CTSE_Info::x_IndexBioseq(CBioseq_Info* info) 
{
    //    x_RegisterRemovedIds(bioseq,info);
    CBioObjectId uniq_id;
    _ASSERT(info->GetBioObjectId().GetType() == CBioObjectId::eUnSet);
    ITERATE ( CBioseq_Info::TId, it, info->GetId() ) {
        if ( it->IsGi() ) {
            uniq_id = CBioObjectId(*it);
            break;
        }
    }
    if ( !info->GetId().empty() ) {
        x_SetBioseqIds(info);
    }
    if (uniq_id.GetType() == CBioObjectId::eUnSet) {
        if (!info->GetId().empty()) 
            uniq_id = CBioObjectId(*info->GetId().begin());
        else 
            uniq_id = x_RegisterBioObject(*info);
    }
    return uniq_id;
}
CBioObjectId CTSE_Info::x_IndexBioseq_set(CBioseq_set_Info* info) 
{
    CBioObjectId uniq_id;
    _ASSERT(info->GetBioObjectId().GetType() == CBioObjectId::eUnSet);
    if (info->m_Bioseq_set_Id > 0)
        uniq_id = CBioObjectId(CBioObjectId::eSetId, info->m_Bioseq_set_Id);
    else 
        uniq_id = x_RegisterBioObject(*info);
    return uniq_id;
}

CBioObjectId CTSE_Info::x_RegisterBioObject(CTSE_Info_Object& info)
{
    CBioObjectId uniq_id = info.GetBioObjectId();
    if (uniq_id.GetType() == CBioObjectId::eUniqNumber) {
        if (m_BioObjects.find(uniq_id) != m_BioObjects.end())
            return uniq_id;
    }
        
    uniq_id = CBioObjectId(CBioObjectId::eUniqNumber,
                           m_InternalBioObjNumber++);
    m_BioObjects[uniq_id] = &info;
    return uniq_id;
}

void CTSE_Info::x_UnregisterBioObject(CTSE_Info_Object& info)
{
    const CBioObjectId& uid = info.GetBioObjectId();
    if (uid.GetType() == CBioObjectId::eUniqNumber) {
        TBioObjects::iterator i = m_BioObjects.find(uid);
        if (i != m_BioObjects.end()) {
            // i->second->x_SetBioObjectId(CBioObjectId());
            m_BioObjects.erase(i);
        }
    }    
}

CTSE_Info_Object* CTSE_Info::x_FindBioObject(const CBioObjectId& uniq_id) const
{
    switch (uniq_id.GetType()) {
    case CBioObjectId::eUniqNumber:
        {
            TBioObjects::const_iterator i = m_BioObjects.find(uniq_id);
            if (i != m_BioObjects.end())
                return i->second;
        }
        return NULL;
    case CBioObjectId::eSeqId:
        {
            x_GetRecords(uniq_id.GetSeqId(), true);
            CFastMutexGuard guard(m_BioseqsMutex);
            TBioseqs::const_iterator it = m_Bioseqs.find(uniq_id.GetSeqId());
            if ( it != m_Bioseqs.end() ) {
                return it->second;
            }
        }
        return NULL;
    case CBioObjectId::eSetId:
        {
            TBioseq_sets::const_iterator it = 
                m_Bioseq_sets.find(uniq_id.GetSetId());
            if ( it != m_Bioseq_sets.end() ) {
                return it->second;
            }
        }
        return NULL;
    default:
        _ASSERT(0);
    }
    return NULL;
}

/*
void CTSE_Info::x_RegisterRemovedIds(const CConstRef<CBioseq>& bioseq, 
                                     CBioseq_Info* info)
{
    if (m_SetObjectInfo) {
        CTSE_SetObjectInfo::TBioseq_InfoMap::iterator iter =
            m_SetObjectInfo->m_Bioseq_InfoMap.find(bioseq);
        if ( iter != m_SetObjectInfo->m_Bioseq_InfoMap.end() ) {
            vector<CSeq_id_Handle>& ids = iter->second.m_Removed_ids;
            for( vector<CSeq_id_Handle>::const_iterator idit = ids.begin();
                 idit != ids.end(); ++idit) {
                _ASSERT( m_Removed_Bioseqs.find(*idit) == m_Removed_Bioseqs.end());
                m_Removed_Bioseqs.insert(TBioseqs::value_type(*idit, info));
            }
        }
    }
}
*/
CRef<CSeq_annot_SNP_Info>
CTSE_Info::x_GetSNP_Info(const CConstRef<CSeq_annot>& annot)
{
    CRef<CSeq_annot_SNP_Info> ret;
    if ( m_SetObjectInfo ) {
        CTSE_SetObjectInfo::TSeq_annot_InfoMap::iterator iter =
            m_SetObjectInfo->m_Seq_annot_InfoMap.find(annot);
        if ( iter != m_SetObjectInfo->m_Seq_annot_InfoMap.end() ) {
            ret = iter->second.m_SNP_annot_Info;
            m_SetObjectInfo->m_Seq_annot_InfoMap.erase(iter);
        }
    }
    return ret;
}


bool CTSE_Info::HasAnnot(const CAnnotName& name) const
{
    TAnnotLockReadGuard guard(GetAnnotLock());
    return m_NamedAnnotObjs.find(name) != m_NamedAnnotObjs.end();
}


bool CTSE_Info::HasNamedAnnot(const string& name) const
{
    return HasAnnot(CAnnotName(name));
}


bool CTSE_Info::HasUnnamedAnnot(void) const
{
    return HasAnnot(CAnnotName());
}


CConstRef<CSeq_entry> CTSE_Info::GetCompleteTSE(void) const
{
    return GetCompleteSeq_entry();
}


CConstRef<CSeq_entry> CTSE_Info::GetTSECore(void) const
{
    return GetSeq_entryCore();
}


void CTSE_Info::x_DSAttachContents(CDataSource& ds)
{
    _ASSERT(!m_DataSource);

    m_DataSource = &ds;
    TParent::x_DSAttachContents(ds);
    if ( m_Split ) {
        m_Split->x_DSAttach(ds);
    }
    ITERATE ( TBioseqs, it, m_Bioseqs ) {
        ds.x_IndexSeqTSE(it->first, this);
    }
    ds.x_IndexAnnotTSEs(this);
}


void CTSE_Info::x_DSDetachContents(CDataSource& ds)
{
    _ASSERT(m_DataSource == &ds);

    ITERATE ( TBioseqs, it, m_Bioseqs ) {
        ds.x_UnindexSeqTSE(it->first, this);
    }
    ds.x_UnindexAnnotTSEs(this);
    TParent::x_DSDetachContents(ds);
    m_DataSource = 0;
}


void CTSE_Info::x_DSMapObject(CConstRef<TObject> obj, CDataSource& ds)
{
    ds.x_Map(obj, this);
    TParent::x_DSMapObject(obj, ds);
}


void CTSE_Info::x_DSUnmapObject(CConstRef<TObject> obj, CDataSource& ds)
{
    ds.x_Unmap(obj, this);
    TParent::x_DSUnmapObject(obj, ds);
}


inline
void CTSE_Info::x_IndexSeqTSE(const CSeq_id_Handle& id)
{
    if ( HasDataSource() ) {
        GetDataSource().x_IndexSeqTSE(id, this);
    }
}


inline
void CTSE_Info::x_UnindexSeqTSE(const CSeq_id_Handle& id)
{
    if ( HasDataSource() ) {
        GetDataSource().x_UnindexSeqTSE(id, this);
    }
}


void CTSE_Info::x_IndexAnnotTSE(const CAnnotName& name,
                                const CSeq_id_Handle& id)
{
    if ( !id.IsGi() ) {
        m_AnnotIdsFlags |= fAnnotIds_NonGi;
        if ( id.HaveMatchingHandles() ) {
            m_AnnotIdsFlags |= fAnnotIds_Matching;
        }
    }
    TIdAnnotInfoMap::iterator iter = m_IdAnnotInfoMap.lower_bound(id);
    if ( iter == m_IdAnnotInfoMap.end() || iter->first != id ) {
        iter = m_IdAnnotInfoMap
            .insert(iter, TIdAnnotInfoMap::value_type(id, SIdAnnotInfo()));
        bool orphan = !ContainsMatchingBioseq(id);
        iter->second.m_Orphan = orphan;
        if ( HasDataSource() ) {
            GetDataSource().x_IndexAnnotTSE(id, this, orphan);
        }
    }
    _VERIFY(iter->second.m_Names.insert(name).second);
}


void CTSE_Info::x_UnindexAnnotTSE(const CAnnotName& name,
                                  const CSeq_id_Handle& id)
{
    TIdAnnotInfoMap::iterator iter = m_IdAnnotInfoMap.lower_bound(id);
    if ( iter == m_IdAnnotInfoMap.end() || iter->first != id ) {
        return;
    }
    _VERIFY(iter->second.m_Names.erase(name) == 1);
    if ( iter->second.m_Names.empty() ) {
        bool orphan = iter->second.m_Orphan;
        m_IdAnnotInfoMap.erase(iter);
        if ( HasDataSource() ) {
            GetDataSource().x_UnindexAnnotTSE(id, this, orphan);
        }
    }
}


void CTSE_Info::x_DoUpdate(TNeedUpdateFlags flags)
{
    if ( flags & (fNeedUpdate_core|fNeedUpdate_children_core) ) {
        if ( m_Split ) {
            m_Split->x_UpdateCore();
        }
    }
    TParent::x_DoUpdate(flags);
}


namespace {
    static inline void x_SortUnique(CTSE_Info::TSeqIds& ids)
    {
        sort(ids.begin(), ids.end());
        ids.erase(unique(ids.begin(), ids.end()), ids.end());
    }
}


void CTSE_Info::GetBioseqsIds(TSeqIds& ids) const
{
    
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        ITERATE ( TBioseqs, it, m_Bioseqs ) {
            ids.push_back(it->first);
        }
    }}
    if ( m_Split ) {
        m_Split->GetBioseqsIds(ids);
        // after adding split bioseq Seq-ids the result may contain
        // duplicates and need to be sorted
        x_SortUnique(ids);
    }
}


void CTSE_Info::GetAnnotIds(TSeqIds& ids) const
{
    {{
        TAnnotLockReadGuard guard(GetAnnotLock());
        ITERATE ( TNamedAnnotObjs, it, m_NamedAnnotObjs ) {
            ITERATE ( TAnnotObjs, it2, it->second ) {
                ids.push_back(it2->first);
            }
        }
    }}
    x_SortUnique(ids);
}


bool CTSE_Info::ContainsBioseq(const CSeq_id_Handle& id) const
{
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        if ( m_Bioseqs.find(id) != m_Bioseqs.end() ) {
            return true;
        }
    }}
    if ( m_Split ) {
        return m_Split->ContainsBioseq(id);
    }
    return false;
}


bool CTSE_Info::ContainsMatchingBioseq(const CSeq_id_Handle& id) const
{
    if ( ContainsBioseq(id) ) {
        return true;
    }
    else if ( id.HaveMatchingHandles() ) {
        CSeq_id_Handle::TMatches ids;
        id.GetMatchingHandles(ids);
        ITERATE ( CSeq_id_Handle::TMatches, match_it, ids ) {
            if ( *match_it != id ) {
                if ( ContainsBioseq(*match_it) ) {
                    return true;
                }
            }
        }
    }
    return false;
}


CConstRef<CBioseq_Info> CTSE_Info::FindBioseq(const CSeq_id_Handle& id) const
{
    CConstRef<CBioseq_Info> ret;
    x_GetRecords(id, true);
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        TBioseqs::const_iterator it = m_Bioseqs.find(id);
        if ( it != m_Bioseqs.end() ) {
            ret = it->second;
        }
    }}
    return ret;
}


CConstRef<CBioseq_Info>
CTSE_Info::FindMatchingBioseq(const CSeq_id_Handle& id) const
{
    return GetSeqMatch(id).m_Bioseq;
}


SSeqMatch_TSE CTSE_Info::GetSeqMatch(const CSeq_id_Handle& id) const
{
    SSeqMatch_TSE ret;
    ret.m_Bioseq = FindBioseq(id);
    if ( ret.m_Bioseq ) {
        ret.m_Seq_id = id;
    }
    else if ( id.HaveMatchingHandles() ) {
        CSeq_id_Handle::TMatches ids;
        id.GetMatchingHandles(ids);
        ITERATE ( CSeq_id_Handle::TMatches, match_it, ids ) {
            if ( *match_it != id ) {
                ret.m_Bioseq = FindBioseq(*match_it);
                if ( ret.m_Bioseq ) {
                    ret.m_Seq_id = *match_it;
                    break;
                }
            }
        }
    }
    return ret;
}


void CTSE_Info::x_GetRecords(const CSeq_id_Handle& id, bool bioseq) const
{
    if ( m_Split ) {
        m_Split->x_GetRecords(id, bioseq);
    }
}


void CTSE_Info::x_LoadChunk(TChunkId chunk_id) const
{
    m_Split->x_LoadChunk(chunk_id);
}


void CTSE_Info::x_LoadChunks(const TChunkIds& chunk_ids) const
{
    m_Split->x_LoadChunks(chunk_ids);
}


void CTSE_Info::x_SetBioseqId(const CSeq_id_Handle& id,
                              CBioseq_Info* info)
{
    _ASSERT(info);
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        pair<TBioseqs::iterator, bool> ins =
            m_Bioseqs.insert(TBioseqs::value_type(id, info));
        if ( !ins.second ) {
            // No duplicate bioseqs in the same TSE
            NCBI_THROW_FMT(CObjMgrException, eAddDataError,
                           "duplicate Bioseq id " << id << " present in" <<
                           "\n  seq1: " << ins.first->second->IdString() <<
                           "\n  seq2: " << info->IdString());
        }
    }}
    // register this TSE in data source as containing the sequence
    x_IndexSeqTSE(id);
}


void CTSE_Info::x_SetBioseqIds(CBioseq_Info* info)
{
    _ASSERT(info);
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        ITERATE ( CBioseq_Info::TId, it, info->GetId() ) {
            pair<TBioseqs::iterator, bool> ins =
                m_Bioseqs.insert(TBioseqs::value_type(*it, info));
            if ( !ins.second ) {
                // No duplicate bioseqs in the same TSE
                NCBI_THROW(CObjMgrException, eAddDataError,
                           "duplicate Bioseq id "+it->AsString()+" present in"+
                           "\n  seq1: " + ins.first->second->IdString()+
                           "\n  seq2: " + info->IdString());
            }
        }
    }}
    // register this TSE in data source as containing the sequence
    if ( HasDataSource() ) {
        GetDataSource().x_IndexSeqTSE(info->GetId(), this);
    }
}


void CTSE_Info::x_ResetBioseqId(const CSeq_id_Handle& id,
                                CBioseq_Info* info)
{
    {{
        CFastMutexGuard guard(m_BioseqsMutex);
        TBioseqs::iterator iter = m_Bioseqs.lower_bound(id);
        if ( iter == m_Bioseqs.end() || iter->first != id ) {
            return;
        }
        _ASSERT(iter->second == info);
        m_Bioseqs.erase(iter);

        if (m_Split) {
            iter = m_Removed_Bioseqs.find(id);
            if (iter == m_Removed_Bioseqs.end())
                m_Removed_Bioseqs.insert(TBioseqs::value_type(id, info));
        }

    }}
    x_UnindexSeqTSE(id);
}


void CTSE_Info::x_SetBioseq_setId(int key,
                                  CBioseq_set_Info* info)
{
    pair<TBioseq_sets::iterator, bool> ins =
        m_Bioseq_sets.insert(TBioseq_sets::value_type(key, info));
    if ( ins.second ) {
        // everything is fine
    }
    else {
        // No duplicate bioseqs in the same TSE
        NCBI_THROW(CObjMgrException, eAddDataError,
                   " duplicate Bioseq_set id '"+NStr::IntToString(key));
    }
}


void CTSE_Info::x_ResetBioseq_setId(int key,
                                    CBioseq_set_Info* info)
{
    TBioseq_sets::iterator iter = m_Bioseq_sets.lower_bound(key);
    if ( iter != m_Bioseq_sets.end() && iter->first == key ) {
        _ASSERT(iter->second == info);
        m_Bioseq_sets.erase(iter);
        if (m_Split) {
            iter = m_Removed_Bioseq_sets.find(key);
            if (iter == m_Removed_Bioseq_sets.end())
                m_Removed_Bioseq_sets.insert(TBioseq_sets::value_type(key, info));
        }
    }    
}


void CTSE_Info::x_SetDirtyAnnotIndexNoParent(void)
{
    if ( HasDataSource() ) {
        GetDataSource().x_SetDirtyAnnotIndex(*this);
    }
}


void CTSE_Info::x_ResetDirtyAnnotIndexNoParent(void)
{
    if ( HasDataSource() ) {
        GetDataSource().x_ResetDirtyAnnotIndex(*this);
    }
}


void CTSE_Info::UpdateAnnotIndex(const CSeq_id_Handle& id) const
{
    x_GetRecords(id, false);
    const_cast<CTSE_Info*>(this)->UpdateAnnotIndex();
}


void CTSE_Info::UpdateAnnotIndex(void) const
{
    const_cast<CTSE_Info*>(this)->UpdateAnnotIndex();
}


void CTSE_Info::UpdateAnnotIndex(const CTSE_Info_Object& object) const
{
    const_cast<CTSE_Info*>(this)->
        UpdateAnnotIndex(const_cast<CTSE_Info_Object&>(object));
}


void CTSE_Info::UpdateAnnotIndex(void)
{
    UpdateAnnotIndex(*this);
}


void CTSE_Info::UpdateAnnotIndex(CTSE_Info_Object& object)
{
    _ASSERT(&object.GetTSE_Info() == this);
    if ( object.x_DirtyAnnotIndex() ) {
        CDataSource::TAnnotLockWriteGuard guard(eEmptyGuard);
        if (HasDataSource())
            guard.Guard(GetDataSource());
        TAnnotLockWriteGuard guard2(GetAnnotLock());
        object.x_UpdateAnnotIndex(*this);
        _ASSERT(!object.x_DirtyAnnotIndex());
    }
}

/*
void CTSE_Info::UpdateAnnotIndex(CTSE_Chunk_Info& chunk)
{
    CDataSource::TAnnotLockWriteGuard guard(GetDataSource());
    TAnnotLockWriteGuard guard2(GetAnnotLock());
    chunk.x_UpdateAnnotIndex(*this);
}
*/

void CTSE_Info::x_UpdateAnnotIndexContents(CTSE_Info& tse)
{
    _ASSERT(this == &tse);
    TParent::x_UpdateAnnotIndexContents(tse);
}


CTSE_Info::TAnnotObjs& CTSE_Info::x_SetAnnotObjs(const CAnnotName& name)
{
    TNamedAnnotObjs::iterator iter = m_NamedAnnotObjs.lower_bound(name);
    if ( iter == m_NamedAnnotObjs.end() || iter->first != name ) {
        typedef TNamedAnnotObjs::value_type value_type;
        iter = m_NamedAnnotObjs.insert(iter, value_type(name, TAnnotObjs()));
    }
    return iter->second;
}


void CTSE_Info::x_RemoveAnnotObjs(const CAnnotName& name)
{
    m_NamedAnnotObjs.erase(name);
}


const CTSE_Info::TAnnotObjs*
CTSE_Info::x_GetAnnotObjs(const CAnnotName& name) const
{
    TNamedAnnotObjs::const_iterator iter = m_NamedAnnotObjs.lower_bound(name);
    if ( iter == m_NamedAnnotObjs.end() || iter->first != name ) {
        return 0;
    }
    return &iter->second;
}


const CTSE_Info::TAnnotObjs*
CTSE_Info::x_GetUnnamedAnnotObjs(void) const
{
    TNamedAnnotObjs::const_iterator iter = m_NamedAnnotObjs.begin();
    if ( iter == m_NamedAnnotObjs.end() || iter->first.IsNamed() ) {
        return 0;
    }
    return &iter->second;
}


SIdAnnotObjs& CTSE_Info::x_SetIdObjects(TAnnotObjs& objs,
                                        const CAnnotName& name,
                                        const CSeq_id_Handle& id)
{
    // repeat for more generic types of selector
    TAnnotObjs::iterator it = objs.lower_bound(id);
    if ( it == objs.end() || it->first != id ) {
        // new id
        it = objs.insert(it, TAnnotObjs::value_type(id, SIdAnnotObjs()));
        x_IndexAnnotTSE(name, id);
    }
    _ASSERT(it != objs.end() && it->first == id);
    return it->second;
}


SIdAnnotObjs& CTSE_Info::x_SetIdObjects(const CAnnotName& name,
                                        const CSeq_id_Handle& id)
{
    return x_SetIdObjects(x_SetAnnotObjs(name), name, id);
}


const SIdAnnotObjs* CTSE_Info::x_GetIdObjects(const TAnnotObjs& objs,
                                              const CSeq_id_Handle& idh) const
{
    TAnnotObjs::const_iterator it = objs.lower_bound(idh);
    if ( it == objs.end() || it->first != idh ) {
        return 0;
    }
    return &it->second;
}


const SIdAnnotObjs* CTSE_Info::x_GetIdObjects(const CAnnotName& name,
                                              const CSeq_id_Handle& idh) const
{
    const TAnnotObjs* objs = x_GetAnnotObjs(name);
    if ( !objs ) {
        return 0;
    }
    return x_GetIdObjects(*objs, idh);
}


const SIdAnnotObjs*
CTSE_Info::x_GetUnnamedIdObjects(const CSeq_id_Handle& idh) const
{
    const TAnnotObjs* objs = x_GetUnnamedAnnotObjs();
    if ( !objs ) {
        return 0;
    }
    return x_GetIdObjects(*objs, idh);
}


bool CTSE_Info::x_HasIdObjects(const CSeq_id_Handle& idh) const
{
    // tse annot index should be locked by TAnnotLockReadGuard
    ITERATE ( TNamedAnnotObjs, it, m_NamedAnnotObjs ) {
        if ( x_GetIdObjects(it->second, idh) ) {
            return true;
        }
    }
    return false;
}


inline
void CTSE_Info::x_MapAnnotObject(TRangeMap& rangeMap,
                                 const SAnnotObject_Key& key,
                                 const SAnnotObject_Index& index)
{
    //_ASSERT(index.m_AnnotObject_Info == key.m_AnnotObject_Info);
    rangeMap.insert(TRangeMap::value_type(key.m_Range, index));
}


inline
bool CTSE_Info::x_UnmapAnnotObject(TRangeMap& rangeMap,
                                   const CAnnotObject_Info& info,
                                   const SAnnotObject_Key& key)
{
    for ( TRangeMap::iterator it = rangeMap.find(key.m_Range);
          it && it->first == key.m_Range; ++it ) {
        if ( it->second.m_AnnotObject_Info == &info ) {
            rangeMap.erase(it);
            return rangeMap.empty();
        }
    }
    _ASSERT(0);
    return rangeMap.empty();
}


void CTSE_Info::x_MapAnnotObject(SIdAnnotObjs& objs,
                                 const SAnnotObject_Key& key,
                                 const SAnnotObject_Index& index)
{
    if ( index.m_AnnotObject_Info->IsLocs() ) {
        // Locs may contain multiple indexes
        CAnnotObject_Info::TTypeIndexSet idx_set;
        index.m_AnnotObject_Info->GetLocsTypes(idx_set);
        ITERATE(CAnnotObject_Info::TTypeIndexSet, idx_rg, idx_set) {
            for (size_t idx = idx_rg->first; idx < idx_rg->second; ++idx) {
                x_MapAnnotObject(objs.x_GetRangeMap(idx), key, index);
            }
        }
    }
    else {
        CAnnotType_Index::TIndexRange idx_rg =
            CAnnotType_Index::GetTypeIndex(*index.m_AnnotObject_Info);
        for (size_t idx = idx_rg.first; idx < idx_rg.second; ++idx) {
            x_MapAnnotObject(objs.x_GetRangeMap(idx), key, index);
        }
    }
}


bool CTSE_Info::x_UnmapAnnotObject(SIdAnnotObjs& objs,
                                   const CAnnotObject_Info& info,
                                   const SAnnotObject_Key& key)
{
    CAnnotType_Index::TIndexRange idx_rg =
        CAnnotType_Index::GetTypeIndex(info);
    for (size_t idx = idx_rg.first; idx < idx_rg.second; ++idx) {
        _ASSERT(idx < objs.x_GetRangeMapCount());
        if ( x_UnmapAnnotObject(objs.x_GetRangeMap(idx), info, key) ) {
            if ( objs.x_CleanRangeMaps() ) {
                return objs.m_SNPSet.empty();
            }
        }
    }
    return false;
}


void CTSE_Info::x_MapAnnotObject(TAnnotObjs& objs,
                                 const CAnnotName& name,
                                 const SAnnotObject_Key& key,
                                 const SAnnotObject_Index& index)
{
    x_MapAnnotObject(x_SetIdObjects(objs, name, key.m_Handle), key, index);
}


bool CTSE_Info::x_UnmapAnnotObject(TAnnotObjs& objs,
                                   const CAnnotName& name,
                                   const CAnnotObject_Info& info,
                                   const SAnnotObject_Key& key)
{
    TAnnotObjs::iterator it = objs.find(key.m_Handle);
    if ( it != objs.end() && x_UnmapAnnotObject(it->second, info, key) ) {
        x_UnindexAnnotTSE(name, key.m_Handle);
        objs.erase(it);
        return objs.empty();
    }
    return false;
}


void CTSE_Info::x_MapSNP_Table(const CAnnotName& name,
                               const CSeq_id_Handle& key,
                               const CSeq_annot_SNP_Info& snp_info)
{
    SIdAnnotObjs& objs = x_SetIdObjects(name, key);
    objs.m_SNPSet.push_back(ConstRef(&snp_info));
}


void CTSE_Info::x_UnmapSNP_Table(const CAnnotName& name,
                                 const CSeq_id_Handle& key,
                                 const CSeq_annot_SNP_Info& snp_info)
{
    SIdAnnotObjs& objs = x_SetIdObjects(name, key);
    TSNPSet::iterator iter = find(objs.m_SNPSet.begin(),
                                  objs.m_SNPSet.end(),
                                  ConstRef(&snp_info));
    if ( iter != objs.m_SNPSet.end() ) {
        objs.m_SNPSet.erase(iter);
    }
}


void CTSE_Info::x_MapAnnotObject(const CAnnotName& name,
                                 const SAnnotObject_Key& key,
                                 const SAnnotObject_Index& index)
{
    x_MapAnnotObject(x_SetAnnotObjs(name), name, key, index);
}


void CTSE_Info::x_UnmapAnnotObject(const CAnnotName& name,
                                   const CAnnotObject_Info& info,
                                   const SAnnotObject_Key& key)
{
    TAnnotObjs& index = x_SetAnnotObjs(name);

    x_UnmapAnnotObject(index, name, info, key);

    if ( index.empty() ) {
        x_RemoveAnnotObjs(name);
    }
}


void CTSE_Info::x_UnmapAnnotObjects(const SAnnotObjectsIndex& infos)
{
    if ( !infos.IsIndexed() ) {
        return;
    }
    const CAnnotName& name = infos.GetName();
    TAnnotObjs& index = x_SetAnnotObjs(name);

    ITERATE ( SAnnotObjectsIndex::TObjectInfos, it, infos.GetInfos() ) {
        if ( it->HasSingleKey() ) {
            x_UnmapAnnotObject(index, name, *it, it->GetKey());
        }
        else {
            for ( size_t i = it->GetKeysBegin(); i < it->GetKeysEnd(); ++i ) {
                x_UnmapAnnotObject(index, name, *it, infos.GetKey(i));
            }
        }
    }

    if ( index.empty() ) {
        x_RemoveAnnotObjs(name);
    }
}


CBioseq_set_Info& CTSE_Info::x_GetBioseq_set(int id)
{
    TBioseq_sets::iterator iter;
    if (m_Split) {
        iter = m_Removed_Bioseq_sets.find(id);
        if ( iter != m_Removed_Bioseq_sets.end() )
            return *iter->second;
    }
        
    iter = m_Bioseq_sets.find(id);
    if ( iter == m_Bioseq_sets.end() ) {
        NCBI_THROW(CObjMgrException, eRegisterError,
                   "cannot find Bioseq-set by local id");
    }
    return *iter->second;
}


CBioseq_Info& CTSE_Info::x_GetBioseq(const CSeq_id_Handle& id)
{
    CFastMutexGuard guard(m_BioseqsMutex);
    TBioseqs::iterator iter;
    if (m_Split) {
        iter = m_Removed_Bioseqs.find(id);
        if ( iter != m_Removed_Bioseqs.end() ) 
            return *iter->second;
    }

    iter = m_Bioseqs.find(id);
    if ( iter == m_Bioseqs.end() ) {
        NCBI_THROW(CObjMgrException, eRegisterError,
                   "cannot find Bioseq by Seq-id "+id.AsString());
    }
    return *iter->second;
}


CTSE_Split_Info& CTSE_Info::GetSplitInfo(void)
{
    if ( !m_Split ) {
        m_Split = new CTSE_Split_Info(GetBlobId(), GetBlobVersion());
        CRef<ITSE_Assigner> listener(new CTSE_Default_Assigner);
        m_Split->x_TSEAttach(*this, listener); 
    }
    return *m_Split;
}


bool CTSE_Info::x_NeedsDelayedMainChunk(void) const
{
    return m_Split && m_Split->x_NeedsDelayedMainChunk();
}


void CTSE_Info::x_AddFeaturesById(TAnnotObjects& objects,
                                  const SFeatIdIndex& index,
                                  TFeatIdInt id,
                                  bool xref) const
{
    if ( !index.m_Chunks.empty() ) {
        x_LoadChunks(index.m_Chunks);
    }
    if ( !index.m_IndexInt ) {
        return;
    }
    const SFeatIdIndex::TIndexInt& index2 = *index.m_IndexInt;
    for ( SFeatIdIndex::TIndexInt::const_iterator iter2 = index2.find(id);
          iter2 != index2.end() && iter2->first == id; ++iter2 ) {
        const SFeatIdInfo& info = iter2->second;
        if ( info.m_Xref == xref ) {
            if ( info.m_IsChunk ) {
                x_LoadChunk(info.m_ChunkId);
            }
            else {
                objects.push_back(info.m_Info);
            }
        }
    }
}


void CTSE_Info::x_AddFeaturesById(TAnnotObjects& objects,
                                  CSeqFeatData::ESubtype subtype,
                                  TFeatIdInt id,
                                  bool xref) const
{
    TFeatIdIndex::const_iterator iter = m_FeatIdIndex.find(subtype);
    if ( iter == m_FeatIdIndex.end() ) {
        return;
    }
    x_AddFeaturesById(objects, iter->second, id, xref);
}


void CTSE_Info::x_AddAllFeaturesById(TAnnotObjects& objects,
                                     TFeatIdInt id,
                                     bool xref) const
{
    //LOG_POST_X(1, this << ": ""x_AddAllFeaturesWithId: " << id);
    ITERATE ( TFeatIdIndex, iter, m_FeatIdIndex ) {
        x_AddFeaturesById(objects, iter->second, id, xref);
    }
}


void CTSE_Info::x_AddFeaturesById(TAnnotObjects& objects,
                                  const SFeatIdIndex& index,
                                  const TFeatIdStr& id,
                                  bool xref) const
{
    if ( !index.m_Chunks.empty() ) {
        x_LoadChunks(index.m_Chunks);
    }
    if ( !index.m_IndexStr ) {
        return;
    }
    const SFeatIdIndex::TIndexStr& index2 = *index.m_IndexStr;
    for ( SFeatIdIndex::TIndexStr::const_iterator iter2 = index2.find(id);
          iter2 != index2.end() && iter2->first == id; ++iter2 ) {
        const SFeatIdInfo& info = iter2->second;
        if ( info.m_Xref == xref ) {
            if ( info.m_IsChunk ) {
                x_LoadChunk(info.m_ChunkId);
            }
            else {
                objects.push_back(info.m_Info);
            }
        }
    }
}


void CTSE_Info::x_AddFeaturesById(TAnnotObjects& objects,
                                  CSeqFeatData::ESubtype subtype,
                                  const TFeatIdStr& id,
                                  bool xref) const
{
    TFeatIdIndex::const_iterator iter = m_FeatIdIndex.find(subtype);
    if ( iter == m_FeatIdIndex.end() ) {
        return;
    }
    x_AddFeaturesById(objects, iter->second, id, xref);
}


void CTSE_Info::x_AddAllFeaturesById(TAnnotObjects& objects,
                                     const TFeatIdStr& id,
                                     bool xref) const
{
    //LOG_POST_X(1, this << ": ""x_AddAllFeaturesWithId: " << id);
    ITERATE ( TFeatIdIndex, iter, m_FeatIdIndex ) {
        x_AddFeaturesById(objects, iter->second, id, xref);
    }
}


CTSE_Info::SFeatIdIndex::TIndexInt&
CTSE_Info::x_GetFeatIdIndexInt(CSeqFeatData::ESubtype type)
{
    //LOG_POST_X(2, this << ": ""x_MapFeatById: " << id << " " << xref);
    SFeatIdIndex& index = m_FeatIdIndex[type];
    if ( !index.m_IndexInt ) {
        index.m_IndexInt.reset(new SFeatIdIndex::TIndexInt);
    }
    return *index.m_IndexInt;
}


CTSE_Info::SFeatIdIndex::TIndexStr&
CTSE_Info::x_GetFeatIdIndexStr(CSeqFeatData::ESubtype type)
{
    //LOG_POST_X(2, this << ": ""x_MapFeatById: " << id << " " << xref);
    SFeatIdIndex& index = m_FeatIdIndex[type];
    if ( !index.m_IndexStr ) {
        index.m_IndexStr.reset(new SFeatIdIndex::TIndexStr);
    }
    return *index.m_IndexStr;
}


void CTSE_Info::x_MapFeatById(TFeatIdInt id,
                              CAnnotObject_Info& info,
                              bool xref)
{
    //LOG_POST_X(2, this << ": ""x_MapFeatById: " << id << " " << xref);
    SFeatIdIndex::TIndexInt& index =
        x_GetFeatIdIndexInt(info.GetFeatSubtype());
    SFeatIdIndex::TIndexInt::value_type value(id, SFeatIdInfo(xref, &info));
    index.insert(value);
}


void CTSE_Info::x_UnmapFeatById(TFeatIdInt id,
                                CAnnotObject_Info& info,
                                bool xref)
{
    //LOG_POST_X(3, this << ": ""x_UnmapFeatById: " << id << " " << xref);
    SFeatIdIndex::TIndexInt& index =
        x_GetFeatIdIndexInt(info.GetFeatSubtype());
    for ( SFeatIdIndex::TIndexInt::iterator iter = index.lower_bound(id);
          iter != index.end() && iter->first == id; ++iter ) {
        if ( iter->second.m_Info == &info && iter->second.m_Xref == xref ) {
            index.erase(iter);
            return;
        }
    }
    _ASSERT("x_UnmapFeatById: not found" && 0);
}


void CTSE_Info::x_MapFeatById(const TFeatIdStr& id,
                              CAnnotObject_Info& info,
                              bool xref)
{
    //LOG_POST_X(2, this << ": ""x_MapFeatById: " << id << " " << xref);
    SFeatIdIndex::TIndexStr& index =
        x_GetFeatIdIndexStr(info.GetFeatSubtype());
    SFeatIdIndex::TIndexStr::value_type value(id, SFeatIdInfo(xref, &info));
    index.insert(value);
}


void CTSE_Info::x_UnmapFeatById(const TFeatIdStr& id,
                                CAnnotObject_Info& info,
                                bool xref)
{
    //LOG_POST_X(3, this << ": ""x_UnmapFeatById: " << id << " " << xref);
    SFeatIdIndex::TIndexStr& index =
        x_GetFeatIdIndexStr(info.GetFeatSubtype());
    for ( SFeatIdIndex::TIndexStr::iterator iter = index.lower_bound(id);
          iter != index.end() && iter->first == id; ++iter ) {
        if ( iter->second.m_Info == &info && iter->second.m_Xref == xref ) {
            index.erase(iter);
            return;
        }
    }
    _ASSERT("x_UnmapFeatById: not found" && 0);
}


void CTSE_Info::x_MapFeatById(const TFeatId& id,
                              CAnnotObject_Info& info,
                              bool xref)
{
    if ( id.IsId() ) {
        x_MapFeatById(id.GetId(), info, xref);
    }
    else {
        x_MapFeatById(id.GetStr(), info, xref);
    }
}


void CTSE_Info::x_UnmapFeatById(const TFeatId& id,
                                CAnnotObject_Info& info,
                                bool xref)
{
    if ( id.IsId() ) {
        x_UnmapFeatById(id.GetId(), info, xref);
    }
    else {
        x_UnmapFeatById(id.GetStr(), info, xref);
    }
}


void CTSE_Info::x_MapFeatByLocus(const string& locus, bool tag,
                                 CAnnotObject_Info& info)
{
    m_LocusIndex.insert(TLocusIndex::value_type(TLocusKey(locus, tag), &info));
}


void CTSE_Info::x_UnmapFeatByLocus(const string& locus, bool tag,
                                   CAnnotObject_Info& info)
{
    for ( TLocusIndex::iterator it =
              m_LocusIndex.lower_bound(TLocusKey(locus, tag));
          it != m_LocusIndex.end() &&
              it->first.first == locus &&
              it->first.second == tag;
          ++it ) {
        if ( it->second == &info ) {
            m_LocusIndex.erase(it);
            return;
        }
    }
}


void CTSE_Info::x_MapChunkByFeatId(TFeatIdInt id,
                                   CSeqFeatData::ESubtype subtype,
                                   TChunkId chunk_id,
                                   bool xref)
{
    SFeatIdIndex::TIndexInt& index = x_GetFeatIdIndexInt(subtype);
    SFeatIdIndex::TIndexInt::value_type value(id, SFeatIdInfo(xref, chunk_id));
    index.insert(value);
}


void CTSE_Info::x_MapChunkByFeatId(TFeatIdInt id,
                                   CSeqFeatData::E_Choice type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    CAnnotType_Index::TIndexRange range =
        CAnnotType_Index::GetFeatTypeRange(type);
    for ( size_t index = range.first; index < range.second; ++index ) {
        CSeqFeatData::ESubtype subtype =
            CAnnotType_Index::GetSubtypeForIndex(index);
        x_MapChunkByFeatId(id, subtype, chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(TFeatIdInt feat_id,
                                   const SAnnotTypeSelector& type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    if ( type.GetFeatSubtype() != CSeqFeatData::eSubtype_any ) {
        x_MapChunkByFeatId(feat_id, type.GetFeatSubtype(), chunk_id, xref);
    }
    else {
        x_MapChunkByFeatId(feat_id, type.GetFeatType(), chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatIdStr& id,
                                   CSeqFeatData::ESubtype subtype,
                                   TChunkId chunk_id,
                                   bool xref)
{
    SFeatIdIndex::TIndexStr& index = x_GetFeatIdIndexStr(subtype);
    SFeatIdIndex::TIndexStr::value_type value(id, SFeatIdInfo(xref, chunk_id));
    index.insert(value);
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatIdStr& id,
                                   CSeqFeatData::E_Choice type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    CAnnotType_Index::TIndexRange range =
        CAnnotType_Index::GetFeatTypeRange(type);
    for ( size_t index = range.first; index < range.second; ++index ) {
        CSeqFeatData::ESubtype subtype =
            CAnnotType_Index::GetSubtypeForIndex(index);
        x_MapChunkByFeatId(id, subtype, chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatIdStr& feat_id,
                                   const SAnnotTypeSelector& type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    if ( type.GetFeatSubtype() != CSeqFeatData::eSubtype_any ) {
        x_MapChunkByFeatId(feat_id, type.GetFeatSubtype(), chunk_id, xref);
    }
    else {
        x_MapChunkByFeatId(feat_id, type.GetFeatType(), chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatId& id,
                                   CSeqFeatData::ESubtype subtype,
                                   TChunkId chunk_id,
                                   bool xref)
{
    if ( id.IsId() ) {
        x_MapChunkByFeatId(id.GetId(), subtype, chunk_id, xref);
    }
    else {
        x_MapChunkByFeatId(id.GetStr(), subtype, chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatId& id,
                                   CSeqFeatData::E_Choice type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    if ( id.IsId() ) {
        x_MapChunkByFeatId(id.GetId(), type, chunk_id, xref);
    }
    else {
        x_MapChunkByFeatId(id.GetStr(), type, chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatId(const TFeatId& feat_id,
                                   const SAnnotTypeSelector& type,
                                   TChunkId chunk_id,
                                   bool xref)
{
    if ( feat_id.IsId() ) {
        x_MapChunkByFeatId(feat_id.GetId(), type, chunk_id, xref);
    }
    else {
        x_MapChunkByFeatId(feat_id.GetStr(), type, chunk_id, xref);
    }
}


void CTSE_Info::x_MapChunkByFeatType(CSeqFeatData::ESubtype subtype,
                                     TChunkId chunk_id)
{
    m_FeatIdIndex[subtype].m_Chunks.push_back(chunk_id);
}


void CTSE_Info::x_MapChunkByFeatType(CSeqFeatData::E_Choice type,
                                     TChunkId chunk_id)
{
    CAnnotType_Index::TIndexRange range =
        CAnnotType_Index::GetFeatTypeRange(type);
    for ( size_t index = range.first; index < range.second; ++index ) {
        CSeqFeatData::ESubtype subtype =
            CAnnotType_Index::GetSubtypeForIndex(index);
        x_MapChunkByFeatType(subtype, chunk_id);
    }
}


void CTSE_Info::x_MapChunkByFeatType(const SAnnotTypeSelector& type,
                                     TChunkId chunk_id)
{
    if ( type.GetFeatSubtype() != CSeqFeatData::eSubtype_any ) {
        x_MapChunkByFeatType(type.GetFeatSubtype(), chunk_id);
    }
    else {
        x_MapChunkByFeatType(type.GetFeatType(), chunk_id);
    }
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::ESubtype subtype,
                             TFeatIdInt id,
                             bool xref) const
{
    UpdateAnnotIndex();
    TAnnotObjects objects;
    if ( subtype == CSeqFeatData::eSubtype_any ) {
        x_AddAllFeaturesById(objects, id, xref);
    }
    else {
        x_AddFeaturesById(objects, subtype, id, xref);
    }
    return objects;
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::E_Choice type,
                             TFeatIdInt id,
                             bool xref) const
{
    UpdateAnnotIndex();
    TAnnotObjects objects;
    if ( type == CSeqFeatData::e_not_set ) {
        x_AddAllFeaturesById(objects, id, xref);
    }
    else {
        CAnnotType_Index::TIndexRange range =
            CAnnotType_Index::GetFeatTypeRange(type);
        for ( size_t index = range.first; index < range.second; ++index ) {
            CSeqFeatData::ESubtype subtype =
                CAnnotType_Index::GetSubtypeForIndex(index);
            x_AddFeaturesById(objects, subtype, id, xref);
        }
    }
    return objects;
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::ESubtype subtype,
                             const TFeatIdStr& id,
                             bool xref) const
{
    UpdateAnnotIndex();
    TAnnotObjects objects;
    if ( subtype == CSeqFeatData::eSubtype_any ) {
        x_AddAllFeaturesById(objects, id, xref);
    }
    else {
        x_AddFeaturesById(objects, subtype, id, xref);
    }
    return objects;
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::E_Choice type,
                             const TFeatIdStr& id,
                             bool xref) const
{
    UpdateAnnotIndex();
    TAnnotObjects objects;
    if ( type == CSeqFeatData::e_not_set ) {
        x_AddAllFeaturesById(objects, id, xref);
    }
    else {
        CAnnotType_Index::TIndexRange range =
            CAnnotType_Index::GetFeatTypeRange(type);
        for ( size_t index = range.first; index < range.second; ++index ) {
            CSeqFeatData::ESubtype subtype =
                CAnnotType_Index::GetSubtypeForIndex(index);
            x_AddFeaturesById(objects, subtype, id, xref);
        }
    }
    return objects;
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::ESubtype subtype,
                             const TFeatId& id,
                             bool xref) const
{
    TAnnotObjects objects;
    if ( id.IsId() ) {
        x_GetFeaturesById(subtype, id.GetId(), xref).swap(objects);
    }
    else {
        x_GetFeaturesById(subtype, id.GetStr(), xref).swap(objects);
    }
    return objects;
}


CTSE_Info::TAnnotObjects
CTSE_Info::x_GetFeaturesById(CSeqFeatData::E_Choice type,
                             const TFeatId& id,
                             bool xref) const
{
    TAnnotObjects objects;
    if ( id.IsId() ) {
        x_GetFeaturesById(type, id.GetId(), xref).swap(objects);
    }
    else {
        x_GetFeaturesById(type, id.GetStr(), xref).swap(objects);
    }
    return objects;
}


CTSE_Info::TAnnotObjects CTSE_Info::x_GetFeaturesByLocus(const string& locus,
                                                         bool tag) const
{
    UpdateAnnotIndex();
    TAnnotObjects objects;
    for ( TLocusIndex::const_iterator it =
              m_LocusIndex.lower_bound(TLocusKey(locus, tag));
          it != m_LocusIndex.end() &&
              it->first.first == locus &&
              it->first.second == tag;
          ++it ) {
        objects.push_back(it->second);
    }
    return objects;
}


CTSE_Info::TSeq_feat_Lock
CTSE_Info::x_FindSeq_feat(const CSeq_id_Handle& loc_id,
                          TSeqPos loc_pos,
                          const CSeq_feat& feat) const
{
    TSeq_feat_Lock ret;
    CSeqFeatData::ESubtype subtype = feat.GetData().GetSubtype();
    size_t index = CAnnotType_Index::GetSubtypeIndex(subtype);
    TRange range(loc_pos, loc_pos);
    ITERATE ( TNamedAnnotObjs, it_n, m_NamedAnnotObjs ) {
        const SIdAnnotObjs* objs = x_GetIdObjects(it_n->second, loc_id);
        if ( !objs ) {
            continue;
        }
        if ( index < objs->x_GetRangeMapCount() &&
             !objs->x_RangeMapIsEmpty(index) ) {
            const TRangeMap& rmap = objs->x_GetRangeMap(index);
            for ( TRangeMap::const_iterator it(rmap.begin(range)); it; ++it ) {
                const CAnnotObject_Info& annot_info =
                    *it->second.m_AnnotObject_Info;
                if ( !annot_info.IsRegular() ) {
                    continue;
                }
                const CSeq_feat* found_feat = annot_info.GetFeatFast();
                if ( found_feat == &feat ) {
                    ret.first.first = &annot_info.GetSeq_annot_Info();
                    ret.second = annot_info.GetAnnotIndex();
                    return ret;
                }
            }
        }
        /*
        if ( subtype == CSeqFeatData::eSubtype_variation &&
             !objs->m_SNPSet.empty() ) {
            
        }
        */
    }
    return ret;
}


CTSE_SetObjectInfo::CTSE_SetObjectInfo(void)
{
}


CTSE_SetObjectInfo::~CTSE_SetObjectInfo(void)
{
}



END_SCOPE(objects)
END_NCBI_SCOPE
