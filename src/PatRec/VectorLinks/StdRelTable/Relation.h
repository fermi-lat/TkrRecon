#ifndef TKRRELATION_H
#define TKRRELATION_H

#include "RelKey.h"
#include <string>
#include <iostream>

/** 
 * @class Relation
 *
 * @brief This class is used to relate pair of objets.
 *
 * The Relation class builts a relationship between two objects. It can be an m to n
 * relation and it stores additional information as a vector of strings
 *
 * This version has been modified from the original by Marco and Riccardo. The big
 * changes are to not require that the related objects derive from DataObject, and, 
 * importantly, that the grouping of Relations is now done with stl lists and maps.
 *
 * @author Marco Frailis
 * @author Riccardo Giannitrapani
 * @author Tracy Usher
 *    
 * $Header: /nfs/slac/g/glast/ground/cvs/Event/Event/RelTable/Relation.h,v 1.3 2002/09/13 04:51:23 frailis Exp $
 */
namespace TkrRecon
{

template <class T1, class T2> class RelTable;

template <class T1, class T2> class Relation
{  
public:
    Relation(T1* obj1, T2* obj2) : m_first(obj1), m_second(obj2) {}
    Relation(T1* obj1, T2* obj2, std::string info);
    Relation(T1* obj1, T2* obj2, std::vector<std::string> infos);

    const T1* getFirst()  const { return m_first.getData(); } 
    T1* getFirst()              { return m_first.getData(); }

    const T2* getSecond() const { return m_second.getData(); } 
    T2* getSecond()             { return m_second.getData(); }    

    /// Add additional information (as a string) to the relation
    void addInfo(std::string inf);
    std::vector<std::string> getInfos() const;

    /// Fill the ASCII output stream
    std::ostream& fillStream( std::ostream& s ) const;

    friend class RelTable<T1,T2>;
    
private:

    // This typedef to define RelListIter to the compiler
    typedef typename RelationList<T1,T2>::RelationListIter RelListIter;

    // Allow RelTable to insert and remove this relation from the list
    void insertInList(RelationList<T1,T2>* list);
    void removeFromList(RelationList<T1,T2>* list);

    // Allow RelTable to insert and remove the first RelKey into its multimap
    void insertFirst(RelKeyMultiMap<T1,T1,T2>* map)  {m_first.insertInMap(map, m_listIter);}
    void removeFirst(RelKeyMultiMap<T1,T1,T2>* map)  {m_first.removeFromMap(map);}
    void setFirst(T1* obj)                           {m_first.setData(obj);}

    // Allow RelTable to insert and remove the second RelKey into its multimap
    void insertSecond(RelKeyMultiMap<T2,T1,T2>* map) {m_second.insertInMap(map, m_listIter);}
    void removeSecond(RelKeyMultiMap<T2,T1,T2>* map) {m_second.removeFromMap(map);}
    void setSecond(T2 *obj)                          {m_second.setData(obj);}
    
    /// Key associated to the first object to be related
    RelKey<T1,T1,T2> m_first;
    /// Key associated to the second object to be related
    RelKey<T2,T1,T2> m_second;
    /// Reference to our position in the Relation list
    RelListIter      m_listIter;
    /// Additional information associated to the relation
    std::vector<std::string> m_infos;
};



template <class T1, class T2> inline Relation<T1,T2>::Relation(T1* obj1, T2* obj2, std::string info): 
  m_first(obj1,this), m_second(obj2,this), m_infos(1,info) {}


template <class T1, class T2> inline Relation<T1,T2>::Relation(T1* obj1, T2* obj2, std::vector<std::string> infos): 
  m_first(obj1,this), m_second(obj2,this), m_infos(infos) {}



template <class T1, class T2> void Relation<T1,T2>::addInfo(std::string inf) 
{
    // Purpose and Method:  This routine add additional information to the relation.
    // Inputs:  inf is the information to be added.
    
  m_infos.push_back(inf);
}

template <class T1, class T2> std::vector<std::string> Relation<T1,T2>::getInfos() const 
{
    // Purpose and Method: This routine get the additional information
    // associated to the relation. 
    // Outputs: a vector of strings

    return m_infos;
}

template <class T1, class T2> void Relation<T1,T2>::insertInList(RelationList<T1,T2>* list)
{
    RelListIter listIter = list->end();
    listIter   = list->insert(listIter, this);
    m_listIter = listIter;

    return;
}     

template <class T1, class T2> void Relation<T1,T2>::removeFromList(RelationList<T1,T2>* list)
{
    list->erase(m_listIter);
    m_listIter = 0;

    return;
}     


template <class T1, class T2> inline std::ostream& Relation<T1,T2>::fillStream( std::ostream& s ) const 
{
    // Fill the ASCII output stream
    s << "      Base class Relation"
      << "\n        First Column: ";
    m_first.toStream(s);
    s << "\n        Second Column: "; 
    m_second.toStream(s);
    s << "\n        Additional Information: ";
  
    std::vector<std::string>::const_iterator i;
    for(i = m_infos.begin(); i != m_infos.end(); i++)
    {
        s << "\n          " << *i;
    }
  
    return s;
}

}
#endif // RELATION_H 

