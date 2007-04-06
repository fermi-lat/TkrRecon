
#ifndef TKRRELKEY_H
#define TKRRELKEY_H

#include <map>
#include <iostream>

/** 
 * @class RelKey
 *
 * @brief This class is used to relate events data.
 *
 * The RelKey class is used by the Relation class to relate events data. It
 * contains four pointers: one points to a particular event data, another
 * points to the next relation containing the same data pointer, one another
 * points to the previous relation and the fourth one points to the first 
 * relation not containing the same data pointer.  This information is useful 
 * to efficiently search for all relations involving a particular event data
 * or to remove a relation.
 * 
 *
 * @author Marco Frailis 
 * @author Riccardo Giannitrapani
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/PatRec/VectorLinks/StdRelTable/RelKey.h,v 1.1 2005/05/26 20:33:06 usher Exp $
 */
namespace TkrRecon
{

// Prototype the Relation class
template <class T1, class T2> class Relation;

// Define here the RelationList class to be a list of pointers to Relations
template <typename T1, typename T2> class RelationList : public std::list<Relation<T1,T2>*> 
{
public:
    typedef typename RelationList<T1,T2>::iterator RelationListIter;
};

// A RelKey will interact with a multimap which will relate the key back to its relation
template <typename T1, typename T2, typename T3> class RelKeyMultiMap : 
          public std::multimap<T1*, typename RelationList<T2,T3>::iterator > {};

// Define the RelKey class here
template <class T1, class T2, class T3> class RelKey 
{    
public:   
    RelKey()       : m_data(0),   m_iterator(0){}
    RelKey(T1* obj): m_data(obj), m_iterator(0){}
    
   ~RelKey() {}
    
    // Provide ability to set and retrieve the "data" or "key"
    void setData(T1* obj) {m_data = obj;}  
    const T1* getData() const { return m_data;}
    T1* getData() {return m_data;}

    /// Fill the ASCII output stream
    void toStream(std::ostream& s) const;

    friend class Relation<T2,T3>;

private:

    // These typedefs will help compiler determine the iterators are a type
    typedef typename RelationList<T2,T3>::iterator      RelationListIter;
    typedef typename RelKeyMultiMap<T1,T2,T3>::iterator RelKeyMultiMapIter;

    // Allow class Relation to "set" and "get" the iterator into the RelKey multimap
    void setMapIter(RelKeyMultiMapIter &relIter) {m_iterator = relIter;}
    RelKeyMultiMapIter getMapIter() const        {return m_iterator;}

    // Allow class Relation to insert and remove this RelKey from the RelKey multimap
    void insertInMap(RelKeyMultiMap<T1,T2,T3>* map, RelationListIter listIter);
    void removeFromMap(RelKeyMultiMap<T1,T2,T3>* map);
    
    /// Pointer to the object to be related
    T1* m_data;
    /// Iterator to its position in the T1 to Relation multimap
    RelKeyMultiMapIter m_iterator;
};

template <class T1, class T2, class T3> 
inline void RelKey<T1,T2,T3>::insertInMap(RelKeyMultiMap<T1,T2,T3>* map, RelationListIter listIter)
{
    RelKeyMultiMapIter mapIter = map->insert(std::pair<T1*, RelationListIter>(m_data, listIter));
    m_iterator = mapIter;
}

template <class T1, class T2, class T3> 
inline void RelKey<T1,T2,T3>::removeFromMap(RelKeyMultiMap<T1,T2,T3>* map)
{
    map->erase(m_iterator);
    m_iterator = 0;
}

template <class T1, class T2, class T3>
inline void RelKey<T1,T2,T3>::toStream(std::ostream& s) const
{
  /// Fill the ASCII output stream
  s  << "\n        Data                    = " << m_data
     << "\n        Iterator                = " << m_iterator;
}

}
#endif // RELKEY_H 
