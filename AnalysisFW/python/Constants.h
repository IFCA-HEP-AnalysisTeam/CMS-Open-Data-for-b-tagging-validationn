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
 const double totalLumi = 2.33; // (/fb) total integrated lumi for 2011 legacy runA
 const double effecLumi = 0.287; // (/pb) effective lumi for --hltpath HLT_Jet60_v*(1-6) of 2011 legacy runA

#endif
