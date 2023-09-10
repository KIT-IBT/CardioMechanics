/*
 **      Name
 **              Material.h
 **
 **
 **      Description
 **              store electrical and mechanical material parameters
 **
 **      History
 **              29.7.93         -fs             creation material.h/pc
 **              22.1.94         -fs             sigma
 **              00.0.98         -br             Cole Cole description
 **              05.1.01         -fs             mechanical properties
 **              20.5.02         -fs             Exceptions
 **              20.5.03         -kc             Tension
 **
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include <kaMachineOS.h>
#include <kaExceptions.h>
#ifdef _STANDARD_C_PLUS_PLUS
# include <complex>
using namespace std;
#else  // ifdef _STANDARD_C_PLUS_PLUS
# include <complex.h>
#endif  // ifdef _STANDARD_C_PLUS_PLUS

static const complex<double> j           = complex<double>(0, 1);
static const double Epsilon0             = 8.854187818e-12;
static const double Mue0                 = 4*M_PI*1E-7;
static const int MaxMechanicalParameters = 24;

enum MaterialConstitutiveLaw {
    NoMaterialLaw, MaterialStVenantKirchhoff, MaterialExponent, MaterialMooneyRivlin, MaterialPoleEquation,
    SpringMassLinear, SpringMassPole, MaterialPolyQ, MaterialExpQ
};

typedef unsigned char MaterialTyp;
class MaterialListe;
class Material;

typedef int CaTensionParameter_Type;
enum CaTensionParameter { maxTensionFiber = 0, maxTensionSheet = 1, maxTensionSheetNormal = 2, nH = 3, Ca50 = 4 };

class vier_werte {
    double tabelle[4];
    
public:
    vier_werte() {
        for (int i = 0; i < 4; i++)
            tabelle[i] = 0;
    }
    
    vier_werte(vier_werte const &a) {
        for (int i = 0; i < 4; i++)
            tabelle[i] = a.tabelle[i];
    }
    
    vier_werte & operator=(vier_werte const from) {
        for (int i = 0; i < 4; i++)
            tabelle[i] = from.tabelle[i];
        return *this;
    }
    
    bool neu(int num, double in) {
        if ((num <= 4) && (num >= 1)) {
            tabelle[num-1] = in;
            return true;
        }
        return false;
    }
    
    double get(int num) const {
        if ((num <= 4) && (num >= 1))
            return tabelle[num-1];
        
        return 0;
    }
    
    friend class Material;
    
    friend ostream & operator<<(ostream &os, const vier_werte out);
    friend istream & operator>>(istream &is, vier_werte &in);
};  // class vier_werte


inline istream & operator>>(istream &is, vier_werte &in) {
    double zwischen;
    
    for (int i = 1; i <= 4; i++) {
        is >> zwischen;
        in.neu(i, zwischen);
    }
    return is;
}

inline ostream & operator<<(ostream &os, const vier_werte out) {
    for (int i = 0; i < 4; i++)
        os << endl << "{ " << i+1 << " } : " << out.tabelle[i];
    os << endl;
    return os;
}

class Extrapolation {
public:
    double ExtrapoliereLinear(double x1, double y1, double x2, double x3, double y3) {
        return (y3-y1)/(x3-x1)*(x2-x1)+y1;
    }
    
    double ExtrapoliereDL(double x1, double y1, double x2, double x3, double y3) {
        return pow(10, ExtrapoliereLinear(log10(x1), log10(y1), log10(x2), log10(x3), log10(y3)));
    }
};


class LFWert {
    double ulFrequenz;
    double dWert;
    
public:
    LFWert(double ulFrequenzArg, double dWertArg) {ulFrequenz = ulFrequenzArg; dWert = dWertArg;}
    
    ~LFWert() {}
    
    double HoleFrequenz() {return ulFrequenz;}
    
    double HoleWert() {return dWert;}
};


class LFWertSubListe : public LFWert {
    LFWertSubListe *pNaechster;
    
public:
    LFWertSubListe(double ulFrequenzArg, double dWertArg) : LFWert(ulFrequenzArg, dWertArg)
    {pNaechster = NULL;}
    
    ~LFWertSubListe() {}
    
    LFWertSubListe *Hinzufuegen(LFWertSubListe *LFWertSubListeArg) {
        LFWertSubListe *pAlt = pNaechster;
        
        pNaechster = LFWertSubListeArg;
        if (LFWertSubListeArg->pNaechster)
            throw kaBaseException("LFWertSubListe::Hinzufuegen - interner Fehler");
        LFWertSubListeArg->pNaechster = pAlt;
        return LFWertSubListeArg;
    }
    
    LFWertSubListe *HoleNaechster() {return pNaechster;}
};


class LFepsilon;
class LFkappa;

class LFWertListe {
    vier_werte delta, tau, alpha;
    double ef, sig;
    
    LFWertSubListe *pErster;
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3) = 0;
    
    LFWertSubListe *Hinzufuegen(LFWert *pLFWertArg, LFWertSubListe *pLFWertSubListeArg) {
        LFWertSubListe *pneu = new LFWertSubListe(pLFWertArg->HoleFrequenz(), pLFWertArg->HoleWert());
        
        assert(pneu);
        
        if (pLFWertSubListeArg) {
            assert(pErster != NULL);
            pLFWertSubListeArg->Hinzufuegen(pneu);
        } else {
            assert(pErster == NULL);
            pErster = pneu;
        }
        return pneu;
    }
    
public:
    LFWertListe() {
        pErster = NULL;
    }
    
    ~LFWertListe() {
        LFWertSubListe *pn = NULL;
        
        for (LFWertSubListe *p = pErster; p != NULL; p = pn) {
            pn = p->HoleNaechster();
            delete p;
        }
    }
    
    virtual double Hole(double ulFrequenzArg) {  // must be virtual for cole-cole
        LFWertSubListe *p, *pa = NULL;
        double dWert;
        
        for (p = pErster; p != NULL; p = p->HoleNaechster()) {
            if (p->HoleFrequenz() >= ulFrequenzArg)
                break;
            pa = p;
        }
        
        if (p && pa)
            dWert =
            Extrapoliere((long unsigned int)pa->HoleFrequenz(), pa->HoleWert(), (long unsigned int)ulFrequenzArg,
                         (long unsigned int)p->HoleFrequenz(), p->HoleWert());
        else if (p)
            dWert = p->HoleWert();
        else if (pa)
            dWert = pa->HoleWert();
        else
            dWert = -1;
        
        return dWert;
    }
    
    LFWertSubListe *HinzufuegenSortiert(LFWert *pLFWertArg) {
        LFWertSubListe *p, *pa = NULL;
        
        for (p = pErster; p != NULL; p = p->HoleNaechster()) {
            if (p->HoleFrequenz() >= pLFWertArg->HoleFrequenz())
                break;
            pa = p;
        }
        return Hinzufuegen(pLFWertArg, pa);
    }
    
    LFWertSubListe *HoleErster() {return pErster;}
    
    // Methods added for cole-cole
    
    complex<double> calc_ecr(double freq) {
        complex<double> ergebnis;
        int lauf;
        complex<double> Disp = 0;
        
        if (freq < 1e-5)
            freq = 1e-5;
        
        for (lauf = 1; lauf <= 4; lauf++) {
            Disp =
            Disp +
            delta.get(lauf)/
            (double(1)+
             pow(complex<double>(j*double(2)*M_PI*freq*tau.get(lauf)),
                 complex<double>((double(1)-alpha.get(lauf)))
                 )
             );
        }
        ergebnis = ef+Disp+(sig/(j*double(2)*M_PI*freq*Epsilon0));
        return ergebnis;
    }
    
    friend class LFepsilon;
    friend class LFkappa;
    friend class Material;
};  // class LFWertListe


class LFmue : public LFWertListe, public Extrapolation {
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3)
    {return ExtrapoliereDL((double)x1, y1, (double)x2, (double)x3, y3);}
    
public:
    LFmue() {}
    
    ~LFmue() {}
};


class LFkappa : public LFWertListe, public Extrapolation {
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3)
    {return ExtrapoliereDL((double)x1, y1, (double)x2, (double)x3, y3);}
    
    bool cole_cole;
    
public:
    double AnisotropyX, AnisotropyY, AnisotropyZ;
    double StretchX, StretchY, StretchZ;
    
    virtual double Hole(double ulFrequenzArg) {  // must be virtual for cole-cole
        LFWertSubListe *p, *pa = NULL;
        double dWert;
        
        if (cole_cole) {
            if (ulFrequenzArg > 1e-3) {
                return -2*M_PI*ulFrequenzArg*Epsilon0*imag(calc_ecr(ulFrequenzArg));
            } else {
                return sig;
            }
        }
        
        for (p = pErster; p; p = p->HoleNaechster()) {
            if (p->HoleFrequenz() >= ulFrequenzArg)
                break;
            pa = p;
        }
        
        if (p && pa) {
            dWert =
            Extrapoliere((long unsigned int)pa->HoleFrequenz(), pa->HoleWert(), (long unsigned int)ulFrequenzArg,
                         (long unsigned int)p->HoleFrequenz(), p->HoleWert());
        } else {
            if (p)
                dWert = p->HoleWert();
            else if (pa)
                dWert = pa->HoleWert();
            else
                dWert = -1;
        }
        
        return dWert;
    }  // Hole
    
    LFkappa() {
        // flag for cole-cole
        cole_cole   = 0; // set to 1 just for the test!!! should be set to 0 !!
        AnisotropyX = AnisotropyY = AnisotropyZ = 1.0;
        StretchX    = StretchY = StretchZ = 0.0;
    }
    
    ~LFkappa() {}
    
    const LFkappa & operator=(LFkappa const &t) {
        AnisotropyX = t.AnisotropyX;
        AnisotropyY = t.AnisotropyY;
        AnisotropyZ = t.AnisotropyZ;
        
        StretchX = t.StretchX;
        StretchY = t.StretchY;
        StretchZ = t.StretchZ;
        return *this;
    }
    
    inline double HoleAnisotropyX() {return AnisotropyX;}
    
    inline double HoleAnisotropyY() {return AnisotropyY;}
    
    inline double HoleAnisotropyZ() {return AnisotropyZ;}
    
    inline double HoleStretchX() {return StretchX;}
    
    inline double HoleStretchY() {return StretchY;}
    
    inline double HoleStretchZ() {return StretchZ;}
    
    friend class Material;
};  // class LFkappa


class LFepsilon  : public LFWertListe, public Extrapolation {
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3)
    {return ExtrapoliereDL((double)x1, y1, (double)x2, (double)x3, y3);}
    
    bool cole_cole;
    
public:
    virtual double Hole(double ulFrequenzArg) {  // must be virtual for cole-cole
        LFWertSubListe *p, *pa = NULL;
        double dWert;
        
        if (cole_cole)
            return real(calc_ecr(ulFrequenzArg));
        
        for (p = LFWertListe::pErster; p != NULL; p = p->HoleNaechster()) {
            if (p->HoleFrequenz() >= ulFrequenzArg)
                break;
            pa = p;
        }
        
        if (p && pa)
            dWert =
            Extrapoliere((long unsigned int)pa->HoleFrequenz(), pa->HoleWert(), (long unsigned int)ulFrequenzArg,
                         (long unsigned int)p->HoleFrequenz(), p->HoleWert());
        else if (p)
            dWert = p->HoleWert();
        else if (pa)
            dWert = pa->HoleWert();
        else
            dWert = -1;
        
        return dWert;
    } // Hole
    
    LFepsilon() {
        // flag for cole-cole
        cole_cole = 0;  // set to 1 just for the test!!! should be set to 0 !!
    }
    
    ~LFepsilon() {}
    
    friend class Material;
};  // class LFepsilon


class LFAttCoeff  : public LFWertListe, public Extrapolation {
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3)
    {return ExtrapoliereDL((double)x1, y1, (double)x2, (double)x3, y3);}
    
public:
    LFAttCoeff() {}
    
    ~LFAttCoeff() {}
};


class LFAbsCoeff  : public LFWertListe, public Extrapolation {
    virtual double Extrapoliere(unsigned long x1, double y1, unsigned long x2, unsigned long x3, double y3)
    {return ExtrapoliereDL((double)x1, y1, (double)x2, (double)x3, y3);}
    
public:
    LFAbsCoeff() {}
    
    ~LFAbsCoeff() {}
};


class LFEModulMaterialParameters {
public:
    MaterialConstitutiveLaw type;
    
    double param[MaxMechanicalParameters];
    
    LFEModulMaterialParameters() {
        type = NoMaterialLaw;
        for (int i = 0; i < MaxMechanicalParameters; i++) param[i] = .0;
    }
    
    ~LFEModulMaterialParameters() {}
};


class CaTension {
    double CaTP[5];
    
public:
    inline CaTension() {
#if KAMATERIALDEBUG
        cerr<<"CaTension::CaTension()"<<endl;
#endif  // if KAMATERIALDEBUG
        for (CaTensionParameter_Type i = maxTensionFiber; i <= Ca50; i++)
            CaTP[i] = 0.;
    }
    
    inline CaTension(double mTF, double mTS, double mTSN, double nHill, double K) {
#if KAMATERIALDEBUG
        cerr<<"CaTension::CaTension(double, double, double, double)"<<endl;
#endif  // if KAMATERIALDEBUG
        CaTP[maxTensionFiber]       = mTF;
        CaTP[maxTensionSheet]       = mTS;
        CaTP[maxTensionSheetNormal] = mTSN;
        CaTP[nH]                    = nHill;
        CaTP[Ca50]                  = K;
    }
    
    ~CaTension()
    {}
    
    inline void Set(CaTensionParameter_Type pp, double value) {
        CaTP[pp] = value;
    }
    
    inline double Get(CaTensionParameter_Type pp) const {
        return CaTP[pp];
    }
    
    inline const CaTension & operator=(CaTension pct) {
        for (CaTensionParameter_Type i = maxTensionFiber; i <= Ca50; i++)
            CaTP[i] = pct.CaTP[i];
        return *this;
    }
};  // class CaTension


class Material : public LFmue, public LFkappa, public LFepsilon, public LFAttCoeff, public LFAbsCoeff,
public LFEModulMaterialParameters, public CaTension {
    char *pName;
    MaterialTyp cKennung;
    double Weight;
    
public:
    Material() {
#if KAMATERIALDEBUG
        printf("Material::Material\n");
#endif  // if KAMATERIALDEBUG
        pName = NULL;
    }
    
    Material(char *name, MaterialTyp kennungval) {
#if KAMATERIALDEBUG
        printf("Material::Material\n");
#endif  // if KAMATERIALDEBUG
        pName    = strdup(name);
        cKennung = kennungval;
        Weight   = 0;
    }
    
    ~Material() {
        free(pName);
    }
    
    const Material & operator=(Material const &t) {
        if (pName)
            free(pName);
        pName                                             = (t.pName ? strdup(t.pName) : NULL);
        cKennung                                          = t.cKennung;
        Weight                                            = t.Weight;
        *dynamic_cast<LFkappa *>(this)                    = *dynamic_cast<const LFkappa *>(&t);
        *dynamic_cast<LFEModulMaterialParameters *>(this) = *dynamic_cast<const LFEModulMaterialParameters *>(&t);
        
        return *this;
    }
    
    inline MaterialTyp HoleKennung() {return cKennung;}
    
    inline char *HoleName() {return pName;}
    
    inline double HoleWeight() {return Weight;}
    
    inline void Read(FILE *fp, MaterialListe *pLM);
};  // class Material


class MaterialSubListe : public Material {
    MaterialSubListe *pNaechster;
    
public:
    MaterialSubListe() {pNaechster = NULL;}
    
    virtual ~MaterialSubListe() {
        assert(pNaechster == NULL);
    }
    
    MaterialSubListe *Hinzufuegen(MaterialSubListe *pMaterialSubListeArg) {
        MaterialSubListe *pAlt = pNaechster;
        
        pNaechster = pMaterialSubListeArg;
        assert(!pMaterialSubListeArg->pNaechster);
        pMaterialSubListeArg->pNaechster = pAlt;
        return pMaterialSubListeArg;
    }
    
    MaterialSubListe *HoleNaechste() {return pNaechster;}
    
    void Aufloesen() {pNaechster = NULL;}
    
    const MaterialSubListe & operator=(Material const &t) {
        Material::operator=(t);
        return *this;
    }
};  // class MaterialSubListe


class MaterialListe {
    MaterialSubListe *pErste;
    
    void Init() {pErste = NULL;}
    
public:
    MaterialListe() {
        Init();
    }
    
    MaterialListe(const char *pDateiName) {
#if KAMATERIALDEBUG
        cerr << "MaterialListe::MaterialListe(" << pDateiName << ")";
#endif  // if KAMATERIALDEBUG
        Init();
        FILE *fp = fopen(pDateiName, "r");
        
        if (!fp)
            throw kaBaseException("MaterialListe::MaterialListe - Can't open file \"%s\"!", pDateiName);
        int c;
        while ((c = fgetc(fp)) != EOF)
            switch (c) {
                case ' ':
                case '\n':
                case '\r':
                case '\t':
                    break;
                    
                case '[': {
                    char name[256];
                    char kennung[256];
                    int  dummy;
                    int  rc = fscanf(fp, "%s %s %d", name, kennung, &dummy);
                    if (rc != 3)
                        throw kaBaseException("MaterialListe::MaterialListe - Error in File %s %s", name, kennung);
                    char kennungval;
                    
                    /*if (strlen(kennung)==1) {
                     kennungval=*kennung;
                     if (kennungval=='L')
                     kennungval=' ';
                     }
                     else*/
                    kennungval = atoi(kennung);
                    
                    Material M(name, kennungval);
                    Material *pM = Hinzufuegen(&M);
                    
                    c = fgetc(fp);
                    if (c != ']')
                        throw kaBaseException("MaterialListe::MaterialListe - ] missing");
                    
                    pM->Read(fp, this);
                    break;
                }
                default: {
                    char buf[80];
                    sprintf(buf, "MaterialListe::MaterialListe - Internal error '%s' char '%c'",
                            pDateiName, c);
                    throw kaBaseException(buf);
                    break;
                }
            }
        fclose(fp);
    }
    
    ~MaterialListe() {
        MaterialSubListe *pn = NULL;
        
        for (MaterialSubListe *p = pErste; p != NULL; p = pn) {
            pn = p->HoleNaechste();
            p->Aufloesen();
            delete p;
        }
    }
    
    Material *Hinzufuegen(Material *pMaterialArg) {
        MaterialSubListe *pneu = new MaterialSubListe();
        
        assert(pneu);
        *pneu = *pMaterialArg;
        if (pErste)
            pErste->Hinzufuegen(pneu);
        else
            pErste = pneu;
        return pneu;
    }
    
    Material *Suchen(char *pName) {
#if KAMATERIALDEBUG
        fprintf(stderr, "MaterialListe::Suchen pName '%c'\n", pName);
#endif  // if KAMATERIALDEBUG
        MaterialSubListe *p;
        
        for (p = pErste; p != NULL; p = p->HoleNaechste())
            if (!strcmp(pName, p->HoleName()))
                break;
        return p;
    }
    
    Material *Suchen(MaterialTyp cKennung) {
#if KAMATERIALDEBUG
        fprintf(stderr, "MaterialListe::Suchen cKennung '%c' %d\n", cKennung, cKennung);
#endif  // if KAMATERIALDEBUG
        MaterialSubListe *p;
        
        for (p = pErste; p != NULL; p = p->HoleNaechste())
            if (cKennung == p->HoleKennung())
                break;
        
        return p;
    }
    
    int Anzahl() {
        int Anzahl = 0;
        
        for (MaterialSubListe *p = pErste; p != NULL; p = p->HoleNaechste())
            Anzahl++;
        return Anzahl;
    }
};  // class MaterialListe

inline void Material::Read(FILE *fp, MaterialListe *pLM) {
    assert(fp);
    
    double f;
    int c;
    char Dimension[40];
    double v;
    
    while ((c = fgetc(fp)) != EOF) {
        if (c == '#') {  // found a reference ??
            unsigned int matnr;
            int rc = fscanf(fp, " %u", &matnr);
            assert(rc == 1);
            Material *pM = pLM->Suchen((MaterialTyp)matnr);
            if (!pM)
                throw kaBaseException("Material::Read - Invalid reference to %d", matnr);
            {
                LFepsilon *pLM = pM;
                for (LFWertSubListe *pLFWSL = pLM->HoleErster();
                     pLFWSL;
                     pLFWSL = pLFWSL->HoleNaechster()) {
                    LFepsilon::HinzufuegenSortiert(pLFWSL);
                }
            }
            {
                LFmue *pLM = pM;
                for (LFWertSubListe *pLFWSL = pLM->HoleErster();
                     pLFWSL;
                     pLFWSL = pLFWSL->HoleNaechster()) {
                    LFmue::HinzufuegenSortiert(pLFWSL);
                }
            }
            {
                LFkappa *pLM = pM;
                for (LFWertSubListe *pLFWSL = pLM->HoleErster();
                     pLFWSL;
                     pLFWSL = pLFWSL->HoleNaechster()) {
                    LFkappa::HinzufuegenSortiert(pLFWSL);
                }
            }
            {
                LFAttCoeff *pLM = pM;
                for (LFWertSubListe *pLFWSL = pLM->HoleErster();
                     pLFWSL;
                     pLFWSL = pLFWSL->HoleNaechster()) {
                    LFAttCoeff::HinzufuegenSortiert(pLFWSL);
                }
            }
            {
                LFAbsCoeff *pLM = pM;
                for (LFWertSubListe *pLFWSL = pLM->HoleErster();
                     pLFWSL;
                     pLFWSL = pLFWSL->HoleNaechster()) {
                    LFAbsCoeff::HinzufuegenSortiert(pLFWSL);
                }
            }
            
            
            for (CaTensionParameter_Type i = maxTensionFiber; i <= Ca50; i++)
                CaTension::Set(i, pM->Get(i));
            
            
            Weight                                            = pM->Weight;
            *dynamic_cast<LFkappa *>(this)                    = *dynamic_cast<LFkappa *>(pM);
            *dynamic_cast<LFEModulMaterialParameters *>(this) = *dynamic_cast<LFEModulMaterialParameters *>(pM);
            
            // Zuweisung der zusaetzlichen Variablen
            this->LFkappa::cole_cole   = pM->LFkappa::cole_cole;
            this->LFepsilon::cole_cole = pM->LFepsilon::cole_cole;
            
            this->LFkappa::delta   = pM->LFkappa::delta;
            this->LFepsilon::delta = pM->LFepsilon::delta;
            this->LFkappa::tau     = pM->LFkappa::tau;
            this->LFepsilon::tau   = pM->LFepsilon::tau;
            this->LFkappa::alpha   = pM->LFkappa::alpha;
            this->LFepsilon::alpha = pM->LFepsilon::alpha;
            this->LFkappa::ef      = pM->LFkappa::ef;
            this->LFepsilon::ef    = pM->LFepsilon::ef;
            this->LFkappa::sig     = pM->LFkappa::sig;
            this->LFepsilon::sig   = pM->LFepsilon::sig;
        } else if ((c != ' ') && (c != '\t') && (c != '\n') && (c != '\r') ) {
            ungetc(c, fp);
            if (c == '[')
                break;
            int rc = fscanf(fp, " %s", Dimension);
            if (rc != 1)
                throw kaBaseException("Material::Read - Error in file");
            
            
            if (!strcmp(Dimension, "TABLE")) {
                // nothing to do everything's alright
            } else if (!strcmp(Dimension, "COLE-COLE")) {
                // I found the signal for cole-cole, set the flag
                
                LFkappa::cole_cole   = 1;
                LFepsilon::cole_cole = 1;
                
                // now read the values for cole-cole
                
                double zwischen;
                int i;
                
                // read the deltas
                
                {
                    fscanf(fp, " %*s");  // skip initial
                }
                
                for (i = 1; i <= 4; i++) {
                    fscanf(fp, " %lE", &zwischen);
                    LFkappa::delta.neu(i, zwischen);  // have to set it for both
                    LFepsilon::delta.neu(i, zwischen);
                }
                
                // read the taus
                
                {
                    fscanf(fp, " %*s");  // skip initial
                }
                
                for (i = 1; i <= 4; i++) {
                    fscanf(fp, " %lE", &zwischen);
                    LFkappa::tau.neu(i, zwischen);  // have to set it for both
                    LFepsilon::tau.neu(i, zwischen);
                }
                
                // read the alphas
                
                {
                    fscanf(fp, " %*s");  // skip initial
                }
                
                for (i = 1; i <= 4; i++) {
                    fscanf(fp, " %lE", &zwischen);
                    LFkappa::alpha.neu(i, zwischen);  // have to set it for both
                    LFepsilon::alpha.neu(i, zwischen);
                }
                
                // read ef
                
                {
                    fscanf(fp, " %*s");  // skip initial
                }
                
                {
                    fscanf(fp, " %lE", &zwischen);
                    LFkappa::ef   = zwischen; // have to set it for both
                    LFepsilon::ef = zwischen;
                }
                
                // read sig
                
                {
                    fscanf(fp, " %*s");  // skip initial
                }
                
                {
                    fscanf(fp, " %lE", &zwischen);
                    LFkappa::sig   = zwischen; // have to set it for both
                    LFepsilon::sig = zwischen;
                }
            } else if (!strcmp(Dimension, "CaTension")) {
                double t;
                for (CaTensionParameter_Type i = maxTensionFiber; i <= Ca50; i++) {
                    int rc = fscanf(fp, " %lE", &t);
                    if (rc == 0)
                        throw kaBaseException("Material::Read - CaTension Parameter less than 5.");
                    CaTension::Set(i, t);
                }
            } else if (!strcmp(Dimension, "Weight")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (Weight)");
                Weight = v;
            } else if (!strcmp(Dimension, "AnisotropyX")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (AnisotropyX)");
                AnisotropyX = v;
            } else if (!strcmp(Dimension, "AnisotropyY")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (AnisotropyY)");
                AnisotropyY = v;
            } else if (!strcmp(Dimension, "AnisotropyZ")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (AnisotropyZ)");
                AnisotropyZ = v;
            } else if (!strcmp(Dimension, "StretchX")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (StretchX)");
                StretchX = v;
            } else if (!strcmp(Dimension, "StretchY")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (StretchY)");
                StretchY = v;
            } else if (!strcmp(Dimension, "StretchZ")) {
                int rc = fscanf(fp, "%lE", &v);
                if (rc != 1)
                    throw kaBaseException("Material::Read - Error in material definition (StretchZ)");
                StretchZ = v;
            } else if (!strcmp(Dimension, "EpsilonR")) {
                int rc = fscanf(fp, "%lE %lE", &f, &v);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in EpsilonR definition");
                LFWert LFWert(f, v *Epsilon0);
                LFepsilon::HinzufuegenSortiert(&LFWert);
            } else if (!strcmp(Dimension, "MueR")) {
                int rc = fscanf(fp, "%lE %lE", &f, &v);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in MueR definition");
                LFWert LFWert(f, v *Mue0);
                LFmue::HinzufuegenSortiert(&LFWert);
            } else if (!strcmp(Dimension, "Kappa")) {
                int rc = fscanf(fp, "%lE %lE", &f, &v);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in Kappa definition");
                LFWert LFWert(f, v);
                LFkappa::HinzufuegenSortiert(&LFWert);
            } else if (!strcmp(Dimension, "AbsCoeff")) {
                int rc = fscanf(fp, "%lE %lE", &f, &v);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in AbsCoeff definition");
                LFWert LFWert(f, v);
                LFAbsCoeff::HinzufuegenSortiert(&LFWert);
            } else if (!strcmp(Dimension, "AttCoeff")) {
                int rc = fscanf(fp, "%lE %lE", &f, &v);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in AttCoeff definition");
                LFWert LFWert(f, v);
                LFAttCoeff::HinzufuegenSortiert(&LFWert);
            } else if (!strcmp(Dimension, "MaterialStVenantKirchhoff")) {
                int rc = fscanf(fp, "%lE %lE", &LFEModulMaterialParameters::param[0], &LFEModulMaterialParameters::param[1]);
                if (rc != 2)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialStVenantKirchhoff;
            } else if (!strcmp(Dimension, "MaterialMooneyRivlin")) {
                int rc = fscanf(fp, "%lE %lE %lE", &LFEModulMaterialParameters::param[0], &LFEModulMaterialParameters::param[1],
                                &LFEModulMaterialParameters::param[2]);
                if (rc != 3)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialMooneyRivlin;
            } else if (!strcmp(Dimension, "MaterialPoleEquation")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0],
                                &LFEModulMaterialParameters::param[1], &LFEModulMaterialParameters::param[2],
                                &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5],
                                &LFEModulMaterialParameters::param[6],
                                &LFEModulMaterialParameters::param[7], &LFEModulMaterialParameters::param[8],
                                &LFEModulMaterialParameters::param[9],
                                &LFEModulMaterialParameters::param[10], &LFEModulMaterialParameters::param[11],
                                &LFEModulMaterialParameters::param[12],
                                &LFEModulMaterialParameters::param[13], &LFEModulMaterialParameters::param[14],
                                &LFEModulMaterialParameters::param[15],
                                &LFEModulMaterialParameters::param[16], &LFEModulMaterialParameters::param[17],
                                &LFEModulMaterialParameters::param[18]
                                );
                if (rc != 19)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialPoleEquation;
            } else if (!strcmp(Dimension, "MaterialExponent")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0], &LFEModulMaterialParameters::param[1],
                                &LFEModulMaterialParameters::param[2], &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5]);
                if (rc != 6)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialExponent;
            } else if (!strcmp(Dimension, "MaterialPolyQ")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0], &LFEModulMaterialParameters::param[1],
                                &LFEModulMaterialParameters::param[2], &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5],
                                &LFEModulMaterialParameters::param[6], &LFEModulMaterialParameters::param[7],
                                &LFEModulMaterialParameters::param[8], &LFEModulMaterialParameters::param[9],
                                &LFEModulMaterialParameters::param[10], &LFEModulMaterialParameters::param[11],
                                &LFEModulMaterialParameters::param[12], &LFEModulMaterialParameters::param[13],
                                &LFEModulMaterialParameters::param[14]);
                if (rc != 15)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialPolyQ;
            } else if (!strcmp(Dimension, "MaterialExpQ")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0], &LFEModulMaterialParameters::param[1],
                                &LFEModulMaterialParameters::param[2], &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5],
                                &LFEModulMaterialParameters::param[6], &LFEModulMaterialParameters::param[7],
                                &LFEModulMaterialParameters::param[8], &LFEModulMaterialParameters::param[9],
                                &LFEModulMaterialParameters::param[10], &LFEModulMaterialParameters::param[11]);
                if (rc != 12)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = MaterialExpQ;
            } else if (!strcmp(Dimension, "SpringMassLinear")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0],
                                &LFEModulMaterialParameters::param[1], &LFEModulMaterialParameters::param[2],
                                &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5],
                                &LFEModulMaterialParameters::param[6],
                                &LFEModulMaterialParameters::param[7], &LFEModulMaterialParameters::param[8],
                                &LFEModulMaterialParameters::param[9]
                                );
                if (rc != 10)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = SpringMassLinear;
            } else if (!strcmp(Dimension, "SpringMassPole")) {
                int rc = fscanf(fp, "%lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE %lE",
                                &LFEModulMaterialParameters::param[0],
                                &LFEModulMaterialParameters::param[1], &LFEModulMaterialParameters::param[2],
                                &LFEModulMaterialParameters::param[3],
                                &LFEModulMaterialParameters::param[4], &LFEModulMaterialParameters::param[5],
                                &LFEModulMaterialParameters::param[6],
                                &LFEModulMaterialParameters::param[7], &LFEModulMaterialParameters::param[8],
                                &LFEModulMaterialParameters::param[9],
                                &LFEModulMaterialParameters::param[10], &LFEModulMaterialParameters::param[11],
                                &LFEModulMaterialParameters::param[12]
                                );
                if (rc != 13)
                    throw kaBaseException("Material::Read - Error in %s definition", Dimension);
                LFEModulMaterialParameters::type = SpringMassPole;
            } else {
                throw kaBaseException("Material::Read Material definition %s unknown", Dimension);
            }
        }
    }
    
#if KAMATERIALDEBUG
    cerr<<"   End of reading material data."<<endl;
#endif  // if KAMATERIALDEBUG
}  // Material::Read

#endif  // ifndef MATERIAL_H
