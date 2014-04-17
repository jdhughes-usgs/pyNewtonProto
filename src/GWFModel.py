#import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import cg, bicgstab


import ReadStructuredxml as sxml

class GWFModel:
    def __init__(self, fxml, structuredModel=True):
        self.structuredModel = structuredModel
        if self.structuredModel:
            self.xmlinput = sxml.StructuredXML(fxml)
            xml_static = self.xmlinput.readStructuredXMLStatic()
        else:
            raise Exception('Unstructured model xml data not supported yet')
        #settings
        self.hclose = xml_static['hclose']
        self.rclose = xml_static['rclose']
        self.inneriterations = xml_static['inneriterations']
        self.outeriterations = xml_static['outeriterations']
        self.newtonraphson = xml_static['newtonraphson']
        #dimensions
        if self.structuredModel:
            self.nlay = xml_static['nlay']
            self.nrow = xml_static['nrow']
            self.ncol = xml_static['ncol']
            self.nrc = self.nrow * self.ncol
            self.neq = self.nlay * self.nrow * self.ncol
            self.dx = xml_static['dx']
            self.dy = xml_static['dy']
            #create connectivity
            self.__createStructuredConnectivity()
            #cell type
            self.celltype = self.__structured3DToUnstructured(xml_static['scelltype'])
            self.head0 = self.__structured3DToUnstructured(xml_static['shead0'])
            #top of model and bottom of each layer
            self.top = self.__structured3DToUnstructured(xml_static['selevations'][0:self.nlay, :, :])
            self.bottom = self.__structured3DToUnstructured(xml_static['selevations'][1:self.nlay+1, :, :])
            #determine cell areas, connection distances, and connection widths for each connection
            self.__createStructuredConnectionDimensions()
            #parameters
            self.kh = self.__structured3DToUnstructured(xml_static['shk'])
        else:
            self.neq = xml_static['neq']
            self.nja = xml_static['nja']
            self.ia = xml_static['ia']
            self.ja = xml_static['ja']
            self.verticalconnection = xml_static['verticalconnection']
            self.cellarea = xml_static['cellarea']
            self.connectionlength = xml_static['connectionlength']
            self.connectionwidth = xml_static['connectionwidth']
            self.celltype = xml_static['celltype']
            self.head0 = xml_static['head0']
            self.top = xml_static['top']
            self.bottom = xml_static['bottom']
        #calculate cell thickness and base Transmissivity for each node
        self.thickness = self.top - self.bottom
        self.Trans = self.kh * self.thickness
        self.cellsaturation = np.ones(self.neq, np.float)
        #create coefficient matrix, rhs, and x
        self.a = np.zeros(self.nja, np.float)
        self.x = np.empty(self.neq, np.float)
        self.rhs = np.zeros(self.neq, np.float)
        #initialize x
        self.x = self.head0[:]

    #update head and boundary condition data for the stress period
    def get_StressPeriodData(self, kper):
        if self.structuredModel:
            data_dict = self.xmlinput.getStressPeriodData(kper)
            if len(data_dict) > 0:
                try:
                    t = data_dict['shead0']
                    self.__set_StructuredHead(data_dict['shead0'])
                except:
                    pass
        else:
            raise Exception('Unstructured model xml data not supported yet')


    #solve
    def solve(self, kper, kstp):
        converged = False
        print 'Solving stress period: {0:5d} time step: {1:5d}'.format(kper+1, kstp+1)
        for outer in xrange(self.outeriterations):
            #assemble matrix
            self.__assemble()
            #solve matrix
            self.acsr = csr_matrix((self.a, self.ja, self.ia), shape=(self.neq, self.neq))
            M = self.get_preconditioner()
            x0 = np.copy(self.x)
            r0 = self.__calculateResidual(x0)
            info = 0
            self.x[:], info = cg(self.acsr, self.rhs, x0=x0, tol=self.rclose, maxiter=self.inneriterations, M=M)
            #self.x[:], info = bicgstab(self.acsr, self.rhs, x0=x0, tol=self.rclose, maxiter=self.inneriterations, M=M)
            if info < 0:
                raise Exception('illegal input or breakdown in linear solver...')
            r1 = self.__calculateResidual(self.x)
            hmax = np.abs(self.x - x0).max()
            if hmax <= self.hclose and abs(r1) <= self.rclose:
                print ' Outer Iterations: {0}'.format(outer+1)
                converged = True
                break
        return converged

    def get_preconditioner(self, **kwargs):
        '''
        Create an ilu preconditioner using spilu.  Return the preconditioner
        as a LinearOperator object.  These are described at:

        http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spilu.html#scipy.sparse.linalg.spilu
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html#scipy.sparse.linalg.LinearOperator

        The idea for this approach was based on information found at:

        http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html

        '''
        drop_tol = None
        fill_factor = None
        drop_rule = None
        permc_spec = None
        diag_pivot_thresh = None
        relax = None
        panel_size = None
        options = None
        if kwargs.has_key('drop_tol'): drop_tol = kwargs['drop_tol']
        if kwargs.has_key('fill_factor'): fill_factor = kwargs['fill_factor']
        if kwargs.has_key('drop_rule'): drop_rule = kwargs['drop_rule']
        if kwargs.has_key('permc_spec'): permc_spec = kwargs['permc_spec']
        if kwargs.has_key('diag_pivot_thresh'): diag_pivot_thresh = kwargs['diag_pivot_thresh']
        if kwargs.has_key('relax'): relax = kwargs['relax']
        if kwargs.has_key('panel_size'): panel_size = kwargs['panel_size']
        if kwargs.has_key('options'): options = kwargs['options']

        #if VERBOSE:
        #    print 'drop_tol: ', drop_tol
        #    print 'fill_factor: ', fill_factor
        #    print 'drop_rule: ', drop_rule
        #    print 'permc_spec: ', permc_spec
        #    print 'diag_pivot_thresh: ', diag_pivot_thresh
        #    print 'relax: ', relax
        #    print 'panel_size: ', panel_size
        #    print 'options: ', options

        from scipy.sparse.linalg import spilu, LinearOperator
        A_ilu = spilu(self.acsr,
            drop_tol=drop_tol,
            fill_factor=fill_factor,
            drop_rule=drop_rule,
            permc_spec=permc_spec,
            diag_pivot_thresh=diag_pivot_thresh,
            relax=relax,
            panel_size=panel_size,
            options=options
            )

        M_x = lambda x: A_ilu.solve(x)
        M = LinearOperator(shape=(self.neq, self.neq), matvec=M_x)
        #M = LinearOperator(shape=(self.nodes, self.nodes), matvec=A_ilu.solve)
        return M

    #assemble matrix
    def __assemble(self):
        self.a.fill(0.0)
        self.rhs.fill(0.0)
        #add conductance
        for node in xrange(self.neq):
            idiag = self.ia[node]
            ia0 = idiag + 1
            ia1 = self.ia[node+1]
            hnode = self.x[node]
            if self.celltype[node] == 0:
#            if self.celltype[node] < 1:
                self.a[idiag] = -1.
                self.rhs[node] = -hnode
                continue
            for jdx in xrange(ia0, ia1):
                nodep = self.ja[jdx]
                hnodep = self.x[nodep]
                if self.celltype[nodep] == 0:
                    continue
                v = self.__calculateConductance(jdx, node, nodep, hnode, hnodep)
                self.a[jdx] = v
                self.a[idiag] -= v
            if self.celltype[node] < 0:
                self.a[idiag] -= 1.e20
                self.rhs[node] -= 1.e20*hnode
        #add boundary conditions
        return

    def __calculateResidual(self, x):
        #assemble matrix
        self.__assemble()
        acsr = csr_matrix((self.a, self.ja, self.ia), shape=(self.neq, self.neq))
        r = acsr.dot(x) - self.rhs
        iloc = np.argmax(np.abs(r))
        rmax = r[iloc]
        return rmax


    def __calculateConductance(self, jdx, node, nodep, h, hp):
        v = 0.0
        if self.verticalconnection[jdx] == 0:
            top = self.top[node]
            bottom = self.bottom[node]
            if h >= top:
                saturation = 1.
            elif h > bottom:
                saturation = (h - bottom) / self.thickness[node]
            else:
                saturation = 0.
            topp = self.top[nodep]
            bottomp = self.bottom[nodep]
            if hp >= topp:
                saturationp = 1.
            elif hp > bottomp:
                saturationp = (hp - bottomp) / self.thickness[nodep]
            else:
                saturationp = 0.
            t = self.Trans[node] * saturation
            tp = self.Trans[nodep] * saturationp
            d = self.connectionlength_n[jdx]
            dp = self.connectionlength_m[jdx]
            w = self.connectionwidth[jdx]
            numer = w * t * tp
            denom = (t * d) + (tp * dp)
            if denom > 0.:
                v = numer / denom
        return v

    #structured data processing functions
    def __createStructuredConnectivity(self):
        self.nodelist = np.zeros((self.nlay, self.nrow, self.ncol), np.int)
        node = 0
        #create node list
        for k in xrange(self.nlay):
            for i in xrange(self.nrow):
                for j in xrange(self.ncol):
                    self.nodelist[k, i, j] = node
                    node += 1
        #node search vectors
        kadd = [-1,0,0,0,0,+1]
        iadd = [0,-1,0,0,+1,+1]
        jadd = [0,0,-1,+1,0,0]
        #create connectivity
        self.ia = []
        self.ja = []
        self.verticalconnection = []
        iapos = 0
        self.ia.append(iapos)
        for k in xrange(self.nlay):
            for i in xrange(self.nrow):
                for j in xrange(self.ncol):
                    #add diagonal
                    node = self.nodelist[k, i, j]
                    self.ja.append(node)
                    self.verticalconnection.append(0)
                    iapos += 1
                    for [kp, ip, jp] in zip(kadd, iadd, jadd):
                        kk = k + kp
                        ii = i + ip
                        jj = j + jp
                        if kk < 0 or kk >= self.nlay:
                            continue
                        if ii < 0 or ii >= self.nrow:
                            continue
                        if jj < 0 or jj >= self.ncol:
                            continue
                        nodep = self.nodelist[kk, ii, jj]
                        self.ja.append(nodep)
                        ivc = 0
                        #upward connection
                        if kk < k:
                            ivc = -1
                        #downward connection
                        elif kk > k:
                            ivc = +1
                        self.verticalconnection.append(ivc)
                        iapos += 1
                    self.ia.append(iapos)
        #convert ia, ja, and verticalconnection to numpy arrays
        self.ia = np.array(self.ia)
        self.ja = np.array(self.ja)
        self.verticalconnection = np.array(self.verticalconnection)
        #determine size of nja and create a matrix
        self.nja = self.ja.shape[0]

    def __createStructuredConnectionDimensions(self):
        self.cellarea = np.zeros(self.neq, np.float)
        self.connectionlength_n = np.zeros(self.nja, np.float)
        self.connectionlength_m = np.zeros(self.nja, np.float)
        self.connectionwidth = np.zeros(self.nja, np.float)
        for node in xrange(self.neq):
            k, i, j = self.__node2lrc(node)
            self.cellarea[node] = self.dx[j] * self.dy[i]
            i0 = self.ia[node] + 1
            i1 = self.ia[node+1]
            for jdx in xrange(i0, i1):
                nodep = self.ja[jdx]
                kk, ii, jj = self.__node2lrc(nodep)
                if kk != k:
                    self.connectionlength_n[jdx] = 0.5 * self.thickness[node]
                    self.connectionlength_m[jdx] = 0.5 * self.thickness[nodep]
                elif ii != i:
                    self.connectionlength_n[jdx] = 0.5 * self.dy[i]
                    self.connectionlength_m[jdx] = 0.5 * self.dy[ii]
                    self.connectionwidth[jdx] = self.dx[j]
                elif jj != j:
                    self.connectionlength_n[jdx] = 0.5 * self.dx[j]
                    self.connectionlength_m[jdx] = 0.5 * self.dx[jj]
                    self.connectionwidth[jdx] = self.dy[i]

    #update state variables
    def __set_StructuredCellType(self, scelltype):
        self.celltype = self.__structured3DToUnstructured(scelltype)

    #update state variables for constant head nodes
    def __set_StructuredHead(self, shead0):
        v = self.__structured3DToUnstructured(shead0)
        for node in xrange(self.neq):
            if self.celltype[node] < 1:
                self.head0[node] = v[node]
                self.x[node] = v[node]

    def __node2lrc(self, node):
        k = node / self.nrc
        nn = node - k * self.nrc
        i = nn / self.nrow
        nn -= i * self.nrow
        j = nn
        return k, i, j

    def __structured2DToUnstructured(self, vs):
        v = np.zeros(self.neq, dtype=vs.dtype)
        for i in xrange(self.nrow):
            for j in xrange(self.ncol):
                node = self.nodelist[0, i, j]
                v[node] = vs[i, j]
        return v

    def __structured3DToUnstructured(self, vs):
        v = np.zeros(self.neq, dtype=vs.dtype)
        for k in xrange(self.nlay):
            for i in xrange(self.nrow):
                for j in xrange(self.ncol):
                    node = self.nodelist[k, i, j]
                    v[node] = vs[k, i, j]
        return v

