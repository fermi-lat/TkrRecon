// $Header: /nfs/slac/g/glast/ground/cvs/TkrRecon/TkrRecon/KalFit.h,v 1.2 2001/02/13 01:50:33 igable Exp $
//----------------------------------------
//
//      Kalman Filter Objects Declarations
//
//               J.A. Hernando, B. Atwood.  Santa Cruz,  7/20/98
//----------------------------------------

#ifndef _KALFIT_H
#define _KALFIT_H 1

#include "TkrRecon/SiClusters.h"


class KalMatrix;
class KalTrack;

//##############
class KalPar
//##############
{
    // Kalman Parameters
public:
    
    // constructors
    KalPar ()
        : position(0), slope(0) 
    {}
    KalPar (double a, double b)
        : position(a), slope(b) 
    {}

    friend KalPar operator* (const KalMatrix&, const KalPar&);
    
    friend double operator*(const KalPar&, const KalPar&);
    friend KalPar operator+(const KalPar&, const KalPar&);
    friend KalPar operator-(const KalPar&, const KalPar&);
    
    inline double getPosition() const {return position;}
    inline double getSlope() const {return slope;}
    
private:
    
    double position;
    double slope;
    
};

double operator*(const KalPar&, const KalPar&);
KalPar operator+(const KalPar&, const KalPar&);	
KalPar operator-(const KalPar&, const KalPar&);

KalPar operator*(const KalMatrix&, const KalPar&);

//##############
class KalMatrix
//###############
{
    // kalman matrices (in fact, normal matrices)
    friend class KalFit;
    
public:
    // constructors
    KalMatrix(){
	a[0][0]=a[0][1]=a[1][0]=a[1][1]=0.;
    }
    // 07/29/99 hsg - added copy constructor
    KalMatrix (const KalMatrix& right) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                a[i][j] = right.a[i][j];
            }
        }
    }
    KalMatrix(double sig2a, double sig2b, double covab){
	a[0][0]=sig2a;
	a[0][1]=a[1][0]=covab;
	a[1][1]=sig2b;
    }
    KalMatrix(double a00,double a01,double a10, double a11){
	a[0][0]=a00;
	a[0][1]=a01;
	a[1][0]=a10;
	a[1][1]=a11;
    }
    
    inline double getsiga() const {return a[0][0];}
    inline double getsigb() const {return a[1][1];}
	inline double getcovab() const {return a[0][1];}
    
    // Manipulators
    static int kfdim;
    
    KalMatrix invert();
    KalMatrix transpose();
    
    friend KalPar operator*(const KalMatrix&, const KalPar&);
    friend KalMatrix operator*(const KalMatrix&, const KalMatrix&);
    friend KalMatrix operator+(const KalMatrix&, const KalMatrix&);
    friend KalMatrix operator-(const KalMatrix&, const KalMatrix&);
    
private:
    
    double a[2][2]; 
    
};

KalMatrix operator*(const KalMatrix&, const KalMatrix&);
KalMatrix operator+(const KalMatrix&, const KalMatrix&);
KalMatrix operator-(const KalMatrix&, const KalMatrix&);

//##############
class KalHit
//##############
{
    
public:
    
    enum TYPE {MEAS, FIT, PRED, SMOOTH, UNKNOWN};
    
    KalHit() 
        : type(UNKNOWN) 
    {}
    KalHit(TYPE t) 
        : type(t)  
    {}
    KalHit(TYPE t, const KalPar& p, const KalMatrix& m)
        :type(t),par(p),cov(m) 
    {}
    KalHit (const KalHit& right)
        : type(right.type), par(right.par), cov(right.cov)
    {}
    
    KalHit changeType(TYPE type);
    
    inline const TYPE& getType() const{return type;}
    inline const KalPar& getPar() const{return par;}
    inline const KalMatrix& getCov() const{return cov;}
    
private:
    
    TYPE type;
    KalPar par;
    KalMatrix cov;
};

//##############
class KalPlane
//##############
{
    
public:
    
    KalPlane() 
        : m_IDHit (-1), m_IDTower(-1), m_IDPlane(-1), m_zplane(0), m_eneplane(0)
    {}
    KalPlane(unsigned id, int kplane, double ene, double z, const KalPar& OrthPar, const KalHit& hit)
        : m_IDPlane(kplane), m_zplane(z), m_eneplane(ene), m_OrthPar(OrthPar) 
    { 
	setIDHit(id);
	setHit(hit); 
    }
    KalPlane(unsigned id, int kplane, double ene, double z, const KalPar& OrthPar)
        : m_IDPlane(kplane), m_zplane(z), m_eneplane(ene), m_OrthPar(OrthPar) 
    {
        setIDHit(id);
    }
    
    // Adding Hits, energy and other variables
    void setHit(const KalHit& hit);
    inline void setEnergy(double const e) {
        m_eneplane = (e < 0.03? 0.03:e);}
    inline void setOrthPar(const KalPar& OrthPar) { m_OrthPar = OrthPar; }
    inline void setIDHit(unsigned id) {
		m_IDHit = id;
		m_IDTower = (int) id/1000000;
    }
    inline void setIDHit(int id, int tower) {
        m_IDHit   = id;
        m_IDTower = tower;
    }
	inline void setIDPlane(int id)  {m_IDPlane = id;}
    inline void setZPlane(double z) {m_zplane = z;}

    // Get Information
    inline unsigned getIDHit()  const{return m_IDHit;}
    inline int getIDTower() const{return m_IDTower;}
    inline int getIDPlane()  const{return m_IDPlane;}
    inline double getZPlane()  const{return m_zplane;}
    inline double getEnergy()  const{return m_eneplane;}
    inline KalPar getOrthPar() const{return m_OrthPar;}
    //	inline KalPlane* getOrthPlane() const{return _mOrthKalPlane;}
    
    // Get Information (compatible with LSQFit)
    KalHit   getHit(KalHit::TYPE type) const;
    Point    getPoint(KalHit::TYPE type) const;
    double   getDeltaChiSq(KalHit::TYPE type) const;
    double   getDeltaChiEne(KalHit::TYPE type) const; 
    double   getSigma(KalHit::TYPE type) const;
    KalMatrix getQmaterial() const {return m_Qmaterial;} 

//    void writeOut(std::ostream& out = std::cout) const;

public:

    // Utility functions for Kalman Filter
    KalHit predicted(KalPlane& kpnext);
    KalHit predicted(KalHit::TYPE typ, double zEnd, int klayerEnd);
    KalHit predicted(KalHit::TYPE typ, int nsteps);
    void setDeltaEne(double ene);
    KalHit filter();
    KalHit smoother(const KalPlane& kplast);
    void clean();   // clean the PRED - FIT - SMOOTH values but not the MEAS
    void clear();   // clean everything
    void removeHit();

    // Static Member functions - Related with the MS introduced by one plane
    //           radlen0=radLen() is the radiation thickness for one plane (theta=0)
    static double theta0ms(double ene, double cosZ, double radlen0);
    static double radLen(int kplane);
    static KalMatrix Q(double ene, double slope, double orth_slope, double radlen0);
 
    
private:
    
    unsigned m_IDHit; // SiCluster Index - code = tower+layer+hit
    int m_IDTower;
    int m_IDPlane;
    //	unsigned m_IDTower;
    
    double m_zplane;
    double m_eneplane;
    
    KalPar m_OrthPar;
    //	KalPlane* _mOrthKalPlane;
    
    //	KalTrack* _mMotherTrack;
    
    KalHit m_hitmeas;
    KalHit m_hitpred;
    KalHit m_hitfit;
    KalHit m_hitsmooth;
    KalMatrix m_Qmaterial;  // covarianve matrix of the material effects 
};
bool operator<(const KalPlane&, const KalPlane&); 
bool operator==(const KalPlane&, const KalPlane&); 


//##############
class KalTrack
//##############
{
    
public:
    
	friend class GlastFit;

    // Constructors
    KalTrack();
    void setIniEnergy(double ene);
    virtual ~KalTrack() {}
    
    // access to the Data - Compatible with LSQFit 
    KalTrack getKalTrack() const {return *this;}
    inline double iniEnergy()              const{return m_energy0;} 

    double positionAtZ(double const z) const;

    inline double position(double deltaZ)  const{return m_x0+deltaZ*m_slopeX;} 
    inline double slope()                  const{return m_slopeX;}
    double errorPosition() const;
    double errorSlope() const;
    double errorSlopeAtVertex() const;
    inline double chiSquare()              const{return m_chisq;}
    inline double chiSquareSmooth()        const{return m_chisqSmooth;}
    inline double KalThetaMS()             const{return m_KalThetaMS;}
    inline double KalEnergy()              const{return m_KalEnergy;}
    inline double scatter()                const{return m_rmsResid;}
    inline int numDataPoints()             const{return(int) kplanelist.size();}
    int numGaps() const;
    inline int numSegmentPoints()          const{return m_numSegmentPoints;}
    inline double chiSquareSegment(double penaltyGap = 0.)  const {return m_chisqSegment + penaltyGap*numGaps();}
    
    double kink(int iplane) const;
    double kinkNorma(int iplane) const;

	// operations
    void clear();
    // void writeOut(std::ostream& out = std::cout) const;

    void draw(gui::DisplayRep& v);

    // LSQ Fit compatible
    double maxResidual(int* index)const;
    Point getHit(unsigned)const;
    unsigned getHitIndex(unsigned)const;
    int compareFits(KalTrack& ktrack); // gives error with const ?!
    
    // Drawing and Printting
    // virtual void printOn(std::ostream &os = std::cout) const;
    void drawTrack(gui::DisplayRep& v);
    void drawTrack(gui::DisplayRep& v, SiCluster::view, KalHit::TYPE);
    void drawChiSq(gui::DisplayRep& v, SiCluster::view, KalHit::TYPE);

    double doFit();
    void filterStep(int iplane);

    double computeChiSqSegment(int nhits, KalHit::TYPE typ = KalHit::SMOOTH);

protected:
    // Do Fit
    void ini();
    KalHit generateFirstFitHit();
    void finish();
    
private:

    // Energy Part
    void eneDetermination();
    
    // segment Part
    int computeNumSegmentPoints(KalHit::TYPE typ = KalHit::SMOOTH);

public:
    
    std::vector<KalPlane> kplanelist;
    
private:
    
    double m_energy0;
    double m_x0;
    double m_slopeX;
    double m_chisq;
    double m_chisqSmooth;
    double m_KalEnergy;
    double m_KalThetaMS;
    double m_rmsResid;
 
    int m_numSegmentPoints;
    double m_chisqSegment;
};

#endif _KALFIT_H

