import sys

#import GWFModel as gwf
import GWFNWTModel as gwf

def main():
    #get name of xml control file
    narg = len(sys.argv)
    iarg = 0
    if narg > 1:
        while iarg < narg-1:
            iarg += 1
            fxml = sys.argv[iarg]

    #--read xml file and return structured model instance
    #model = gwf.GWFModel(fxml)
    model = gwf.GWFNWTModel(fxml)

    number_periods = 2
    number_timesteps = 1

    for kper in xrange(number_periods):
        model.get_StressPeriodData(kper)
        for kstp in xrange(number_timesteps):
            converge = model.solve(kper, kstp)
            if converge:
                break
        for t in model.x:
            print t


if __name__ == '__main__':

    main()
