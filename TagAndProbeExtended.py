from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel


class TagAndProbeExtended(PhysicsModel):

    def __init__(self):
        super(PhysicsModel, self).__init__()

    def _getProcessCategory(self, process):
        ret = None
        if self.DC.isSignal[process] or True:
            for cat in self._categories:
                if cat in process:
                    ret = cat
                    break
        return ret

    def setPhysicsOptions(self, physOptions):
        for po in physOptions[:]:
            if po.startswith("categories="):  # shorthand:  categories=cat1,cat2,cat3
                physOptions.remove(po)
                self._categories = po.replace("categories=", "").split(",")
        super(TagAndProbeExtended, self).setPhysicsOptions(physOptions)

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        pois = []
        for cat in self._categories:
            self.modelBuilder.doVar("SF_%s[1,0.,2]" % cat)
            pois.append('SF_%s' % cat)
        self.modelBuilder.doSet("POI", ','.join(pois))

        exp_pass = {}
        exp_fail = {}
        for b in self.DC.bins:
            for p in self.DC.exp[b].keys():
                cat = self._getProcessCategory(p)
                if cat is None:
                    continue
                if 'pass' in b:
                    exp_pass[cat] = self.DC.exp[b][p]
                elif 'fail' in b:
                    exp_fail[cat] = self.DC.exp[b][p]

        for cat in self._categories:
            self.modelBuilder.factory_('expr::fail_scale_{cat}("max(0.,({pass_exp}+{fail_exp}-({pass_exp}*@0))/{fail_exp})", SF_{cat})'.format(cat=cat, pass_exp=exp_pass[cat], fail_exp=exp_fail[cat]))


    def getYieldScale(self, bin, process):
        "Return the name of a RooAbsReal to scale this yield by or the two special values 1 and 0 (don't scale, and set to zero)"
        if self.DC.isSignal[process] or True :
            cat = self._getProcessCategory(process)
            if cat is None:
                return 1
            if 'pass' in bin:
                return 'SF_%s' % cat
            else:
                return 'fail_scale_%s' % cat
        else:
            return 1


tagAndProbe = TagAndProbeExtended()
