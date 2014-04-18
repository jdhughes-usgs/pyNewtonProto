import sys
import os
import xml.etree.ElementTree as xml
import numpy as np

class StructuredXML:
    def __init__(self, xml_file):
        self.xml_file = xml_file

        try:
            tree = xml.parse(self.xml_file)
        except IOError:
            print '{0} does not exist'.format(os.path.basename(xml_file))
        #--parse xml file
        self.root = tree.getroot()

    def readStructuredXMLStatic(self):
        dimensions = self.root.find('Dimensions')
        self.nlay = int(dimensions.find('nlay').text)
        self.nrow = int(dimensions.find('nrow').text)
        self.ncol = int(dimensions.find('ncol').text)
        #--dx
        clist = self.__extractList(handle=dimensions, tag='dx')
        dx = self.__listTo1DArray(clist, shape=(self.ncol,))
        #--dy
        clist = self.__extractList(handle=dimensions, tag='dy')
        dy = self.__listTo1DArray(clist, shape=(self.nrow,))
        #--properties
        properties = self.root.find('Properties')
        #--cell types (<0 constant, 0 inactive, >0 active)
        child = properties.find('CellType')
        scelltype = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol), dtype=np.int)
        #--top and bottom elevations
        selevations = np.zeros((self.nlay+1, self.nrow, self.ncol), np.float)
        #--top
        clist = self.__extractList(handle=properties, tag='TopElevation')
        selevations[0, :, :] = self.__listTo2DArray(clist, shape=(self.nrow, self.ncol))
        #--bottom
        child = properties.find('BottomElevation')
        selevations[1:, :, :] = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol))
        #--Horizontal conductivity
        child = properties.find('HorizontalConductivity')
        shk = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol))
        #--initial state
        initial = self.root.find('InitialState')
        child = initial.find('Head')
        shead0 = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol))
        #--settings
        settings = self.root.find('Settings')
        outeriterations = int(settings.find('OuterIterations').text)
        inneriterations = int(settings.find('InnerIterations').text)
        hclose = float(settings.find('Hclose').text)
        rclose = float(settings.find('Rclose').text)
        try:
            lt = settings.find('NewtonRaphson').text
        except:
            lt = None
        if lt == 'True':
            newtonraphson = True
        else:
            newtonraphson = False
        try:
            lt = settings.find('ConstantHeadAsGHB').text
        except:
            lt = None
        if lt == 'True':
            chasghb = True
        else:
            chasghb = False
        #--return dictionary of data
        return {'nlay':self.nlay, 'nrow':self.nrow, 'ncol':self.ncol, 'dx':dx, 'dy':dy,
                'scelltype':scelltype, 'shead0':shead0,
                'selevations':selevations, 'shk':shk,
                'outeriterations':outeriterations, 'inneriterations':inneriterations,
                'hclose':hclose, 'rclose':rclose,
                'newtonraphson':newtonraphson,
                'chasghb':chasghb}

    def getStressPeriodData(self, kper):
        if kper == 0:
            self.stressperiod = self.root.find('StressPeriodData')
        data_dict = {}
        onPeriod = self.stressperiod.find('StressPeriod{0}'.format(kper+1))
        if onPeriod is None:
            print 'Using stress period data from last stress period'
        else:
            #head data
            child = onPeriod.find('Head')
            if child is not None:
                data_dict['shead0'] = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol))
            #boundary data
            #recharge data
            child = onPeriod.find('Recharge')
            if child is not None:
                data_dict['srecharge'] = self.__get3DDataFromList(handle=child, shape=(self.nlay, self.nrow, self.ncol))
        return data_dict




    def __get3DDataFromList(self, handle=None, shape=None, dtype=np.float):
        if handle is None:
            raise Exception('Valid handle not provided.')
        if len(shape) != 3:
            raise Exception('Three dimensional shape expected.')
        v = np.zeros(shape, dtype=dtype)
        for k in xrange(shape[0]):
            clist = self.__extractList(handle=handle, tag='Layer{0}'.format(k+1))
            v[k, :, :] = self.__listTo2DArray(clist, shape=(shape[1], shape[2]), dtype=dtype)
        return v


    def __extractList(self, handle=None, tag=None):
        if handle is None:
            raise Exception('extractList: handle not provided')
        if tag is None:
            raise Exception('extractList: tag not provided')
        node = handle.find(tag)
        if node is None:
            cerr = '"{0}" tag not found in xml file'.format(tag)
            raise Exception(cerr)
        return node.text.split()


    def __listTo1DArray(self, clist, shape=(1,), dtype=np.float):
        if len(shape) != 1:
            cerr = 'One-dimensional shape expected. Actual shape {1} provided.'.format(shape)
            raise Exception(cerr)
        v = np.ones(shape, np.float)
        if len(clist) == 1:
            if dtype == np.float:
                v *= float(clist[0])
            elif dtype == np.int:
                v *= int(clist[0])
        else:
            vt = []
            for c in clist:
                if dtype == np.float:
                    vt.append(float(c))
                elif dtype == np.int:
                    vt.append(int(c))
            v = np.array(vt)
            if v.shape[0] != shape[0]:
                cerr = 'list shape {0} exceeds specified shape {1}'.format(v.shape, shape)
                raise Exception(cerr)
        return v


    def __listTo2DArray(self, clist, shape=(1, 1), dtype=np.float):
        if len(shape) != 2:
            cerr = 'Two-dimensional shape expected. Actual shape {1} provided.'.format(shape.shape)
            raise Exception(cerr)
        v = np.ones(shape, np.float)
        if len(clist) == 1:
            if dtype == np.float:
                v *= float(clist[0])
            elif dtype == np.int:
                v *= int(clist[0])
        else:
            vt = []
            for c in clist:
                if dtype == np.float:
                    vt.append(float(c))
                elif dtype == np.int:
                    vt.append(int(c))
            v = np.array(vt)
            v = np.reshape(v, shape)
            if v.shape != shape:
                cerr = 'list shape {0} is not equal to specified shape {1}'.format(v.shape, shape)
                raise Exception(cerr)
        return v


if __name__ == '__main__':
    
    narg = len(sys.argv)
    iarg = 0
    if narg > 1:
        while iarg < narg-1:
            iarg += 1
            fxml = sys.argv[iarg]
    sxml = StructuredXML(fxml)
    smdata = sxml.readStructuredXMLStatic()
    print smdata.keys()
    print '...finished'
