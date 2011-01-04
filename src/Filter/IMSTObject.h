/**
 * @class MinSpanTree
 *
 * @brief This defines an interface to the base objects operated on by the MST and Cluster Analyses
 *
 * @author Tracy Usher
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/src/Filter/MinSpanTree.h,v 1.0 2010/12/16 20:44:46 usher Exp $
 */

#ifndef IMSTObject_h
#define IMSTObject_h

#include <vector>
#include <list>
#include <set>
#include <map>

// Forward declaration
class Point;

// Define an interface class so we can operate on several different kinds of objects
class IMSTObject
{
public:
    // Make sure destructor called...
    virtual ~IMSTObject() {}

    // Our object must be able to return the bilayer it is associated with
    virtual const int    getBiLayer()                     const = 0;
    // And, of course, our object must be able to return its position
    virtual const Point& getPosition()                    const = 0;
    // Objects need to determine the distance to a neighbor
    virtual const double getDistanceTo(const IMSTObject&) const = 0;
};

// Some typedefs which enable one to define the adjacency list for the above
typedef std::pair<const IMSTObject*, const IMSTObject*>    MSTObjectPair;
typedef std::pair<MSTObjectPair, double>                   MSTObjectDistPair;
typedef std::map<MSTObjectPair, double>                    MSTObjectToDistMap;
typedef std::list<MSTObjectDistPair >                      MSTObjectDistPairList;

typedef std::pair<const IMSTObject*, double>               MSTObjectToDistPair;
typedef std::pair<const IMSTObject*, MSTObjectToDistPair > MSTObjectToObjectToDistPair;
typedef std::map<const IMSTObject*, double>                MSTObjectDistMap;
typedef std::map<const IMSTObject*, MSTObjectDistMap >     MSTObjectToObjectDistMap;

// More useful stuff
typedef std::vector<IMSTObject*>                           MSTObjectVec;
typedef std::set<const IMSTObject*>                        MSTObjectSet;

#endif

