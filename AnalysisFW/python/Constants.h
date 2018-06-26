#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1


// Flavour selection
// ------------------------------------------------------------------------------------
  enum {
         allflavour,
         b_quark,
         c_quark,
         lgluon,
         b_gsplitting,
         nflavour,
        };

 const TString sflavour [nflavour+1] = {
         "allflavours",
         "b_quark",
         "c_quark",
         "lgluon",
         "b_gsplitting",
         "data"
        };

// Goblal values
// ------------------------------------------------------------------------------------
 const double integratedLumi = 2.33; // (/fb) 2011 legacy runA

#endif
