
#ifndef _H_TkrEventModel_
#define _H_TkrEventModel_

/* Definition of the event structure in the Temporary Transient Data Store.
 *  
 * Only two levels in the logical path are foreseen at present, 
 *   /event/<namespace>/<leave>  e.g. /Event/MC/McVertices
 * 
 * Convention:
 *  If the <leave> object is a
 *  DataObject    use name of corresponding class
 *  Container     use name of ContainedObject class in plural
 *                or append 'Vec' to the name, e.g. use
 *                McVertices or McVertexVec
 *                
 *
 * @author : adapted from LHCb EventModel
 * @author   I. Gable
 * @author   S. Gillespie
 * @author   T. H.-Kozanecka
 */ 

#include <string>

#if defined(_TkrEventModel_CPP_)
#define  _EXTERN_ 
#else
#define  _EXTERN_ extern
#endif

    namespace TkrEventModel {
        _EXTERN_ std::string   EventHeader;

        namespace MC {
            _EXTERN_ std::string Event;
//            _EXTERN_ std::string McEventStructure;
//            _EXTERN_ std::string McPartToHitTab;
//            _EXTERN_ std::string McClusToLyrHitTab;
//            _EXTERN_ std::string McLyrToHitTab;
//            _EXTERN_ std::string McSiLayerHitCol;
            _EXTERN_ std::string PatHitToLyrHit;
            _EXTERN_ std::string PatCandToMcCand;
            _EXTERN_ std::string McPatCandCol;
        }
    }

#undef _EXTERN_
#endif // GLASTEVENT_EVENTMODEL_H
