#ifndef TKRRELTABLE_H
#define TKRRELTABLE_H

#include "Relation.h"
#include <vector>
#include <algorithm>
#include <iterator>

/** 
* @class RelTable
*
* @brief This class is used to wrap a collection of Relations.
*
* The RelTable class wraps a list (a Gaudi ObjectList) of Relations. It
* lets the user search for all object related to a given one. The search can be
* done with respect to an object of the first or the second field of the
* relations. The user can also modify or delete relations.
* 
*
* @author Marco Frailis
* @author Riccardo Giannitrapani
*   
* $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/StdRelTable/RelTable.h,v 1.1 2005/05/26 20:33:06 usher Exp $
*/
namespace TkrRecon 
{

template <class T1, class T2> class RelTable 
{    
public: 
    RelTable() : m_relations(0), m_firstMMap(0), m_secondMMap(0) {}
  
    /// Initialize the internal pointer to an ObjectList of relations
    void init(); // { m_relations = new std::list< Relation<T1,T2>* >;}

    /**
    * The following method add a Relation to the table if it doesn't contain
    * a relation between the same two objects, otherwise it appends the info
    * vector to the exsisting relation
    * @param rel is a pointer to a relation between two objects
    * @return true if the relation has been added and false if it is a duplicate
    * and has not been added (in this case the user has to delete it)
    */
    bool addRelation(Relation<T1,T2>* rel);
  
    /**
    * This method search for all relations having obj in the first
    * field. 
    * @param obj it's a pointer to the object given by the user
    * @return A vector of pointers to the relations involving the given object.
    */
    std::vector< Relation<T1,T2>* > getRelByFirst(const T1* pobj) const;
  
    /**
    * This method search for all relations having pobj in the second
    * field. 
    * @param pobj it's a pointer to the object given by the user
    * @return A vector of pointers to the relations involving the given object.
    */
    std::vector< Relation<T1,T2>* > getRelBySecond(const T2* pobj) const;
  
    /**
    * This method erase a particular relation from the table (keeping the 
    * integrity).
    * @param rel it's a pointer to the relation to be erased
    */
    void erase(Relation<T1,T2> *rel);
  
    /**
    * This method clears the Relation list and the multimaps
    */
    void clear();
  
    /// This method returns the number of relations in the table
    unsigned long size() const ;
  
    /// Returns the pointer to the collection of relations.
    RelationList<T1,T2>* getAllRelations() const;
  
private:

    typedef typename RelationList<T1,T2>::iterator RelationListIter;

    typedef typename RelKeyMultiMap<T1,T1,T2>::iterator    mapT1RelIter;
    typedef typename std::pair<mapT1RelIter, mapT1RelIter> mapT1RelIterPair;

    typedef typename RelKeyMultiMap<T2,T1,T2>::iterator    mapT2RelIter;
    typedef typename std::pair<mapT2RelIter, mapT2RelIter> mapT2RelIterPair;
  
    /// Pointer to a collection of relations
    RelationList<T1,T2>*      m_relations;

    /// Pointer to T1 multimap
    RelKeyMultiMap<T1,T1,T2>* m_firstMMap;

    /// Pointer to T2 multimap
    RelKeyMultiMap<T2,T1,T2>* m_secondMMap;
};
  

template <class T1, class T2> inline void RelTable<T1,T2>::init()
{
    m_relations  = new RelationList<T1,T2>;
    m_firstMMap  = new RelKeyMultiMap<T1,T1,T2>;
    m_secondMMap = new RelKeyMultiMap<T2,T1,T2>;
}

template <class T1,class T2> bool RelTable<T1,T2>::addRelation(Relation<T1,T2>* rel) 
{
    // Purpose and Method:  This routine add a relation to the table if it doesn't 
    // contain a relation between the same two objects, otherwise it appends the info
    // vector to the exsisting relation
    // Inputs:  rel is a pointer to the relation to be added.
    // Outputs: a boolean value which is true if the realtion has been added to the
    //          table and false it it is a duplicate and thus has not been added.
    //          In the latter case the user has to delete the relation

    // This adds the relation to the list
    rel->insertInList(m_relations);

    // Take care of the maps
    rel->insertFirst(m_firstMMap);
    rel->insertSecond(m_secondMMap);

    return true;
}

template <class T1,class T2>
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelByFirst(const T1* pobj) const 
{
    // Purpose and Method: This routine finds all relations having pobj in the
    // first field.  
    // Inputs: pobj is a pointer to the object to be searched in the first field
    // Outputs: A pointer to a vector of Relation* including pobj
    
    std::vector< Relation<T1,T2>* > rels;
    if (!m_relations->size()) return rels;

    T1* pObjLocal = const_cast<T1*>(pobj);

    mapT1RelIterPair iterPair = m_firstMMap->equal_range(pObjLocal);

    for(mapT1RelIter mapIter = iterPair.first; mapIter != iterPair.second; mapIter++)
    {
        RelationList<T1,T2>::RelationListIter relIter = (*mapIter).second;
        Relation<T1,T2>* relation = *relIter;
        rels.push_back(relation);
    }

    return rels;
} 
  
template <class T1,class T2>
    std::vector< Relation<T1,T2>* > RelTable<T1,T2>::getRelBySecond(const T2* pobj) const 
{
    // Purpose and Method: This routine finds all relations having pobj in the
    // second field.  
    // Inputs: pobj is a pointer to the object to be searched in the second field
    // Outputs: A pointer to a vector of Relation* including pobj
    std::vector< Relation<T1,T2>* > rels;
    if (!m_relations->size()) return rels;

    T2* pObjLocal = const_cast<T2*>(pobj);

    mapT2RelIterPair iterPair = m_secondMMap->equal_range(pObjLocal);

    for(mapT2RelIter mapIter = iterPair.first; mapIter != iterPair.second; mapIter++)
    {
        RelationList<T1,T2>::RelationListIter relIter = (*mapIter).second;
        Relation<T1,T2>* relation = *relIter;
        rels.push_back(relation);
    }

    return rels;
}
  
template <class T1,class T2> void RelTable<T1,T2>::erase(Relation<T1,T2>* rel) 
{
    // Purpose: This method remove the given relation from the table

    rel->removeFirst(m_firstMMap);
    rel->removeSecond(m_secondMMap);
    rel->removeFromList(m_relations);

    delete rel;
}
  
template <class T1,class T2> void RelTable<T1,T2>::clear() 
{
    // Purpose: This method remove the given relation from the table

    m_firstMMap->clear();
    m_secondMMap->clear();

    for(RelationList<T1,T2>::RelationListIter relIter = m_relations->begin(); relIter != m_relations->end(); relIter++)
    {
        delete *relIter;
    }

    m_relations->clear();

    return;
}


template <class T1,class T2>
    inline unsigned long RelTable<T1,T2>::size() const 
{
    // Purpose: This method returns the total number of relations contained in the
    // collection
    return m_relations->size();
}
  
template <class T1,class T2> inline RelationList<T1,T2>* RelTable<T1,T2>::getAllRelations() const 
{
    return m_relations;
}

}
#endif // RELTABLE_H 
