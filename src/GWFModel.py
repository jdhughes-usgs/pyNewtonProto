import math
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
        self.rclosetype = xml_static['rclosetype']
        self.inneriterations = xml_static['inneriterations']
        self.outeriterations = xml_static['outeriterations']
        self.newtonraphson = xml_static['newtonraphson']
        self.numericalderiv = xml_static['numericalderiv']
        self.backtracking = xml_static['backtracking']
        self.bottomflag = xml_static['bottomflag']
        self.headsolution = xml_static['headsolution']
        self.averaging = xml_static['averaging']
        self.upw = xml_static['upw']
        self.quadsfactor = xml_static['quadsfactor']
        self.chasghb = xml_static['chasghb']
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
            #set intial heads that are below the bottom of a node to the bottom of the node
            #idx = self.head0 < self.bottom
            #self.head0[idx] = self.bottom[idx] + 1e-6
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
        #boundary condition flags
        self.rechargebnd = False
        #create coefficient matrix, rhs, and x
        self.a = np.zeros(self.nja, np.float)
        self.x = np.empty(self.neq, np.float)
        self.rhs = np.empty(self.neq, np.float)
        self.r = np.empty(self.neq, np.float)
        #initialize x
        self.x = self.head0[:]

    #update head and boundary condition data for the stress period
    def get_StressPeriodData(self, kper):
        if self.structuredModel:
            data_dict = self.xmlinput.getStressPeriodData(kper)
            if len(data_dict) > 0:
                #heads
                try:
                    t = data_dict['shead0']
                    self.__set_StructuredHead(data_dict['shead0'])
                except:
                    pass
                #boundaries
                #recharge
                try:
                    t = data_dict['srecharge']
                    self.__set_StructuredRecharge(data_dict['srecharge'])
                    self.rechargebnd = True
                except:
                    pass
        else:
            raise Exception('Unstructured model xml data not supported yet')


    #solve
    def solve(self, kper, kstp):
        converged = False
        print 'Solving stress period: {0:5d} time step: {1:5d}'.format(kper+1, kstp+1)
        for outer in xrange(self.outeriterations):
            #--create initial x (x0) from a copy of x
            x0 = np.copy(self.x)
            #--assemble conductance matrix
            self.__assemble()
            #--create sparse matrix for residual calculation and conductance formulation
            self.acsr = csr_matrix((self.a, self.ja, self.ia), shape=(self.neq, self.neq))
            #--save sparse matrix with conductance
            if self.newtonraphson:
                self.ccsr = self.acsr.copy()
            else:
                self.ccsr = self.acsr
            #--calculate initial residual
            #  do not attempt a solution if the initial solution is an order of
            #  magnitude less than rclose
            rmax0 = self.__calculateResidual(x0)
            #if outer == 0 and abs(rmax0) <= 0.1 * self.rclose:
            #    break
            if self.backtracking:
                l2norm0 = np.linalg.norm(self.r)
            if self.newtonraphson:
                self.__assemble(nr=True)
                self.acsr = csr_matrix((self.a, self.ja, self.ia), shape=(self.neq, self.neq))
                b = -self.r.copy()
                aif self.headsolution:
                    t = self.acsr.dot(x0)
                    b += t
                else:
                    self.x.fill(0.0)
            else:
                b = self.rhs.copy()
            #--construct the preconditioner
            #M = self.get_preconditioner(fill_factor=3, drop_tol=1e-4)
            M = self.get_preconditioner(fill_factor=3, drop_tol=1e-4)
            #--solve matrix
            info = 0
            if self.newtonraphson:
                self.x[:], info = bicgstab(self.acsr, b, x0=self.x, tol=self.rclose, maxiter=self.inneriterations, M=M)
            else:
                self.x[:], info = cg(self.acsr, b, x0=self.x, tol=self.rclose, maxiter=self.inneriterations, M=M)
            if info < 0:
                raise Exception('illegal input or breakdown in linear solver...')
            #--add upgrade to x0
            #if self.newtonraphson:
            #    if not self.headsolution:
            #        self.x += x0
            #--calculate updated residual
            rmax1 = self.__calculateResidual(self.x)
            #
            if self.bottomflag:
                self.adjusthead(self.x)
            #--back tracking
            if self.backtracking and rmax1 > self.rclose:
                l2norm1 = np.linalg.norm(self.r)
                if l2norm1 > 0.99 * l2norm0:
                    if self.headsolution:
                        dx = self.x - x0
                    else:
                        dx = self.x
                    lv = 0.99
                    for ibk in xrange(100):
                        self.x = x0 + lv * dx
                        rt = self.__calculateResidual(self.x, reset_ccsr=True)
                        rmax1 = rt
                        l2norm = np.linalg.norm(self.r)
                        if l2norm < 0.90 * l2norm0:
                            break
                        lv *= 0.95
            #--calculate hmax
            hmax = np.abs(self.x - x0).max()
            #--calculate
            if hmax <= self.hclose and abs(rmax1) <= self.rclose:
                print ' Outer Iterations: {0}'.format(outer+1)
                converged = True
                self.__calculateQNodes(self.x)
                #print self.cellQ[4], self.cellQ[-4]
                break
        return converged

    def adjusthead(self, x):
        for idx in xrange(self.neq):
            if x[idx] < self.bottom[idx]:
                x[idx] = self.bottom[idx] + 1.e-6
        return

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
    def __assemble(self, x=None, nr=False):
        self.a.fill(0.0)
        self.rhs.fill(0.0)
        if x is None:
            x = self.x
        #--add conductance
        for node in xrange(self.neq):
            idiag = self.ia[node]
            ia0 = idiag + 1
            ia1 = self.ia[node+1]
            hnode = x[node]
            #--check if this cell is inactive or constant
            if self.celltype[node] == 0:
                self.a[idiag] = -1.
                self.rhs[node] = -hnode
                continue
            elif self.celltype[node] < 0:
                if not self.chasghb:
                    self.a[idiag] -= 1.e20
                    self.rhs[node] -= 1.e20*hnode
                else:
                    self.a[idiag] = -1.
                    self.rhs[node] = -hnode
                    continue
            for jdx in xrange(ia0, ia1):
                nodep = self.ja[jdx]
                hnodep = x[nodep]
                if self.celltype[nodep] == 0:
                    continue
                if nr:
                    dx = self.__get_perturbation(hnode)
                    if self.numericalderiv:
                        dh = (hnodep - hnode)
                        q1 = self.__calculateConductance(jdx, node, nodep, hnode, hnodep) * dh
                        dh = (hnodep - (hnode + dx))
                        q2 = self.__calculateConductance(jdx, node, nodep, hnode, hnodep, dx=dx) * dh
                        v = (q1 - q2) / dx
                    else:
                        c1 = self.__calculateConductance(jdx, node, nodep, hnode, hnodep)
                        dcdx = self.__calculateConductance(jdx, node, nodep, hnode, hnodep, anald=True)
                        dq = c1 * dx - (dcdx * dx * (hnodep - (hnode + dx)))
                        v = dq / dx
                else:
                    v = self.__calculateConductance(jdx, node, nodep, hnode, hnodep)
                if self.celltype[nodep] > 0:
                    self.a[idiag] -= v
                    self.a[jdx] = v
                elif self.celltype[nodep] < 0:
                    self.a[idiag] -= v
                    if self.chasghb:
                        self.rhs[node] -= v*hnodep
                    else:
                        self.a[jdx] = v
        #add boundary conditions
        if self.rechargebnd:
            for node in xrange(self.neq):
                if self.celltype[node] > 0:
                    self.rhs[node] -= self.recharge[node]
        return

    def __get_perturbation(self, v):
        #return 1.0e-7
        return np.sqrt(np.finfo(float).eps)

    def __calculateResidual(self, x, reset_ccsr=False):
        #assemble matrix
        self.__assemble(x=x)
        if reset_ccsr:
            self.ccsr = csr_matrix((self.a, self.ja, self.ia), shape=(self.neq, self.neq))
        self.r = self.ccsr.dot(x) - self.rhs
        if self.rclosetype == 'infinity':
            iloc = np.argmax(np.abs(self.r))
            rmax = self.r[iloc]
        elif self.rclosetype == 'l2norm':
            rmax = np.linalg.norm(self.r)
        return rmax

    def __calculateQNodes(self, x):
        #assemble matrix
        self.__assemble(x=x)
        self.cellQ = np.zeros(self.nja, np.float)
        for node in xrange(self.neq):
            idiag = self.ia[node]
            i0 = idiag + 1
            i1 = self.ia[node+1]
            for jdx in xrange(i0, i1):
                jcol = self.ja[jdx]
                self.cellQ[jdx] = self.a[jdx] * (x[jcol] - x[node])
        return None


    def __calculateConductance(self, jdx, node, nodep, h, hp, dx=0.0, anald=False):
        frac1 = 1.005
        frac2 = 0.995
        v = 0.0
        if self.verticalconnection[jdx] == 0:
            if h >= hp:
                hup = h
                tup = self.thickness[node]
                idxh = 0
            else:
                hup = hp
                tup = self.thickness[nodep]
                idxh = 1
            h += dx
            top = self.top[node]
            bottom = self.bottom[node]
            if h >= top:
                #saturation = 1.
                b = self.thickness[node]
            elif h > bottom:
                #saturation = (h - bottom) / self.thickness[node]
                b = (h - bottom)
            else:
                #saturation = 0.
                b = 0.
            topp = self.top[nodep]
            bottomp = self.bottom[nodep]
            if hp >= topp:
                #saturationp = 1.
                bp = self.bottom[nodep]
            elif hp > bottomp:
                #saturationp = (hp - bottomp) / self.thickness[nodep]
                bp = (hp - bottomp)
            else:
                #saturationp = 0.
                bp = 0.
            if idxh == 0:
                bup = hup - bottom
                if bup > self.thickness[node]:
                    bup = self.thickness[node]
            elif idxh == 1:
                bup = hup - bottomp
                if bup > self.thickness[nodep]:
                    bup = self.thickness[nodep]
            if bup < 0.:
                bup = 0.
            k = self.kh[node]
            kp = self.kh[nodep]
            d = self.connectionlength_n[jdx]
            dp = self.connectionlength_m[jdx]
            w = self.connectionwidth[jdx]

            if self.upw:
                hv = k
                hvp = kp
            else:
                hv = k * b
                hvp = kp * bp
            #harmonic
            if self.averaging == 0:
                v = w * hv * hvp / ((hv * dp) + (hvp * d))
            #logarithmic-mean
            elif self.averaging == 1:
                ratio = hvp / hv
                if ratio > frac1 or ratio < frac2:
                    v = (hvp - hv) / math.log(ratio)
                else:
                    v = 0.5 * (hvp + hv)
                v *= w / (dp + d)
            #arithmetic-mean thickness and logarithmic-mean hydraulic conductivity
            elif self.averaging == 2:
                ratio = kp / k
                if ratio > frac1 or ratio < frac2:
                    v = (kp - k) / math.log(ratio)
                else:
                    v = 0.5 * (kp + k)
                v *= w / (dp + d)

            if self.upw:
                bratio = bup / tup
                if anald:
                    a = 1. / (1 - self.quadsfactor)
                    if bratio <= 0.0:
                        sf = 0.
                    elif bratio > 0. and bratio <= self.quadsfactor:
                        sf = a * bratio / (self.quadsfactor * tup)
                    elif bratio > self.quadsfactor and bratio <= (1. - self.quadsfactor):
                        sf = a / tup
                    elif bratio > (1. - self.quadsfactor) and bratio < 1.:
                        sf = 1. - (-1. * a * (1. - bratio)) / (self.quadsfactor * tup)
                    else:
                        sf = 0.
                else:
                    if bratio <= 0.0:
                        v = 1e-6
                        sf = 1.
                        tup = 1.
                    else:
                        a = 1. / (1 - self.quadsfactor)
                        if bratio <= self.quadsfactor:
                            sf = 0.5 * a * (bratio**2) / self.quadsfactor
                        elif bratio > self.quadsfactor and bratio <= (1. - self.quadsfactor):
                            sf = a * bratio + 0.5 * (1. - a)
                        elif bratio > (1. - self.quadsfactor) and bratio < 1.:
                            sf = 1. - ((0.5 * a * ((1. - bratio)**2)) / self.quadsfactor)
                        else:
                            sf = 1.
                v *= tup * sf
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

    #update state variables for recharge boundary
    def __set_StructuredRecharge(self, srecharge):
        self.recharge = self.__structured3DToUnstructured(srecharge)
        self.recharge *= self.cellarea


    def __node2lrc(self, node):
        k = node / self.nrc
        nn = node - k * self.nrc
        i = nn / self.ncol
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

