/*
****************************************************
				SimMedTK LICENSE
****************************************************

****************************************************
*/

#ifndef SMMEMORYBLOCK_H
#define SMMEMORYBLOCK_H

#include "smCore/smConfig.h"
#include "smCore/smCoreClass.h"
#include "smUtilities/smVec3.h"
#include "smCore/smErrorLog.h"
#include <QHash>

enum smMemReturnType{
	SOFMIS_MEMORY_ALLOCATED,
	SOFMIS_MEMORY_ALREADYALLOCATED,
	SOFMIS_MEMORY_NOTENOUGHMEMORY,
	SOFMIS_MEMORY_MEMORYFOUND,
	SOFMIS_MEMORY_NOMEMORYFOUND,
	SOFMIS_MEMORY_INVALIDPARAMS,
	SOFMIS_MEMORY_INVALIDMEMORY,
	SOFMIS_MEMORY_NOERROR
};

///Memory Block makes easy to allocate and associate particular memory.
///it facilities the memory mamangement
class smMemoryBlock:public smCoreClass{

private:
	smErrorLog *log;
	QHash<QString ,void*> memoryBlocks;

public:
	///constructr needs logger in case
	smMemoryBlock(smErrorLog *log){
		type=SOFMIS_SMMEMORYBLOCK;		
		this->log=log;
	}

	smMemoryBlock(){
		type=SOFMIS_SMMEMORYBLOCK;		
		//this->log=smSDK::getErrorLog();
	}

	///alocate a class and returns p_returnedBlock as allocated memory and 
	///return params are SOFMIS_MEMORY_ALLOCATED or SOFMIS_MEMORY_ALREADYALLOCATED or SOFMIS_MEMORY_INVALIDMEMORY
	template<class T> 
	smMemReturnType allocate(QString &p_memoryBlockName,T**p_returnedBlock){
		*p_returnedBlock=new T();
		if(p_returnedBlock==NULL)
			return	SOFMIS_MEMORY_NOTENOUGHMEMORY;
		
		if(memoryBlocks.contains(p_memoryBlockName)){
			delete [] *p_returnedBlock;
			return SOFMIS_MEMORY_ALREADYALLOCATED;
		}
		else{
			memoryBlocks[p_memoryBlockName]=*p_returnedBlock;
			return	SOFMIS_MEMORY_ALLOCATED;
		}
	}

	///alocate any c;asses and returns p_returnedBlock as allocated memory and 
	///return params are SOFMIS_MEMORY_ALLOCATED or SOFMIS_MEMORY_ALREADYALLOCATED or SOFMIS_MEMORY_INVALIDMEMORY
	template<class T> 
	smMemReturnType allocate(QString &p_memoryBlockName,smInt nbr,T**p_returnedBlock){
		*p_returnedBlock=new T[nbr];
		if(p_returnedBlock==NULL)
			return	SOFMIS_MEMORY_NOTENOUGHMEMORY;
		
		if(memoryBlocks.contains(p_memoryBlockName)){
			delete [] *p_returnedBlock;
			return SOFMIS_MEMORY_ALREADYALLOCATED;
		}
		else{
			memoryBlocks[p_memoryBlockName]=*p_returnedBlock;
			return	SOFMIS_MEMORY_ALLOCATED;
		}
	}

	///alocate any classes and returns p_returnedBlock as allocated memory
	///returns	 SOFMIS_MEMORY_ALLOCATED or SOFMIS_MEMORY_ALREADYALLOCATED or SOFMIS_MEMORY_INVALIDMEMORY
	template <class T>
	smMemReturnType allocate(const QString& p_memoryBlockName, const smInt &p_nbr){
		T *allocatedMem;
		allocatedMem=new T[p_nbr];
		if(allocatedMem==NULL)
			return	SOFMIS_MEMORY_NOTENOUGHMEMORY;
		
		if(memoryBlocks.contains(p_memoryBlockName)){
			delete []allocatedMem;
			return SOFMIS_MEMORY_ALREADYALLOCATED;
		}
		else{
			memoryBlocks[p_memoryBlockName]=allocatedMem;
			return	SOFMIS_MEMORY_ALLOCATED;
		}
	}

	///alocate vectors and returns	SOFMIS_MEMORY_INVALIDPARAMS or SOFMIS_MEMORY_ALLOCATED based on block size given
	virtual smMemReturnType allocate(const QString &p_memoryBlockName, const smInt blockSize,void **p_returnedBlock){
		if(blockSize<=0)
			return SOFMIS_MEMORY_INVALIDPARAMS;
		*p_returnedBlock=new smChar[blockSize];
		memoryBlocks[p_memoryBlockName]=*p_returnedBlock;
		return SOFMIS_MEMORY_ALLOCATED;
	}

	///deletes the block from memeory as well as in has container
	virtual smMemReturnType deleteMemory(QString & p_memoryBlockName){
		void *memoryBlock;
		if(memoryBlocks.contains(p_memoryBlockName)){
			memoryBlock=memoryBlocks[p_memoryBlockName];
			delete []memoryBlock;
			memoryBlocks.remove(p_memoryBlockName);
			return SOFMIS_MEMORY_NOERROR;
		}
		else
			return SOFMIS_MEMORY_NOMEMORYFOUND;
	}

	///gets  memory from the container via given block name. it returns SOFMIS_MEMORY_MEMORYFOUND or SOFMIS_MEMORY_NOMEMORYFOUND. 
	virtual smMemReturnType getBlock(const QString &p_memoryBlockName,void **p_memoryPointer){
		if(memoryBlocks.contains(p_memoryBlockName)){
			*p_memoryPointer=memoryBlocks[p_memoryBlockName];
			return	SOFMIS_MEMORY_MEMORYFOUND;
		}
		else
			return SOFMIS_MEMORY_NOMEMORYFOUND;
	}

	template<class T>
	smMemReturnType localtoOriginalBlock(const QString &p_memoryBlockName, T* dst, const smInt p_nbr){

		T *src;
		if(dst!=NULL){
			src=(T*)memoryBlocks[p_memoryBlockName];
			if(src!=NULL){
		 		memcpy(dst,src,sizeof(T)*p_nbr);
				return	SOFMIS_MEMORY_NOERROR;
			}
			else
				 return SOFMIS_MEMORY_INVALIDMEMORY;
		}
		else
			return SOFMIS_MEMORY_INVALIDPARAMS;
	}

	template<class T>
	smMemReturnType originaltoLocalBlock(const QString &p_memoryBlockName, T *src,
	                                     const smInt p_nbr){

		T *dst;
		if(src!=NULL){
			dst=(T*)memoryBlocks[p_memoryBlockName];
			if(dst!=NULL){
		 		memcpy(dst,src,sizeof(T)*p_nbr);
				return	SOFMIS_MEMORY_NOERROR;
			}else
				return	SOFMIS_MEMORY_INVALIDMEMORY;
		}
		else
			return SOFMIS_MEMORY_INVALIDPARAMS;
	}

	///lists the blocks within the container
	void listofBlocks(){
		QHash<QString ,void*>::iterator iter = memoryBlocks.begin();
		while (iter != memoryBlocks.end() ){
			cout << iter.value() << endl;
			++iter;
		}
	}

	//test the destructor
	~smMemoryBlock(){
		foreach (void*memoryPtr, memoryBlocks)
			delete []memoryPtr;
	}

};

#endif
