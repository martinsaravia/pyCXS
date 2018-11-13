import numpy as np, re, time
import msalg as ma
import os
clc = lambda: os.system('cls')
clc()

class CrossSection(object):
    """ Section File reading class...
    """
    def __init__(self, name):
#        pathfile = os.path.dirname(os.path.abspath(__file__))
        cxsdir = os.getcwd() + '\\' + 'inputs' + '\\'
        cxstime = time.time()
        self.readit(cxsdir + name)
        self.pregen()
#        self.offset(0.0,0.0)
        self.vectex()
        self.writex(name)
        self.layers()
        self.center()
        self.stiff()
        print '-----> CXS Total Elapsed Time was:', time.time() - cxstime, ' seconds'

    def readit(self, name):
        start_time = time.time()
        """ Read input file function """
#        print '- Running input file: ', os.getcwd()
        print '-----> Running CXS for input file:', name
        mateid = 0
        lamiid = 0
        maxlam = 0  # Maximum number of laminas
        secfile = open( (name + '.sinp'), 'r')
        linea = secfile.readline()
        elemlist = []
        nodelist = []
        matelist = []
        lamilist = []
        laminames = []
        matenames = []
        while  linea[0:4] != '*END':
            #========================== READ NODE DATA ===================
            if linea[0:5] == '*NODE':
                linea = secfile.readline()
                while linea[0] != '*':
                    linealist = map(float, linea.split(','))
                    linealist[0] = int(linealist[0])
                    nodelist.append(linealist)
                    linea = secfile.readline()

            #========================== READ ELEMENT DATA ===================
            elif linea[0:8] == '*ELEMENT':
#                Store Element Information
                linea = linea.replace(" ", "")
                eset = re.search('(?<=ELSET=)\w+', linea).group()
                lset = re.search('(?<=LAMINATE=)\w+', linea).group()
                etyp = re.search('(?<=TYPE=)\w+', linea).group()
                linea = secfile.readline() # Read New Line
                while linea[0] != '*':
                    linealist = map(int, linea.split(','))
                    linealist.append(eset)
                    linealist.append(lset)
                    linealist.append(etyp)
                    elemlist.append(linealist)
                    linea = secfile.readline()

            elif linea[0:2] == '**':
                linea = secfile.readline()

            elif linea[0:9] == '*MATERIAL':
                linea = linea.replace(" ", "")
                mateid = mateid + 1
                matename = re.search('(?<=NAME=)\w+', linea).group()
                linea = secfile.readline()
                if linea[0:8] == '*ELASTIC':
                    linea = linea.replace(" ", "")
                    matetype = re.search('(?<=TYPE=)\w+', linea).group()
                else:
                    print '- ERROR: Unknown Material', linea,
                    print 'Sorry, exiting...'
                    break
                linea = secfile.readline()
                matecons = map(float, linea.split(','))
                linealist = [matename, matetype, matecons]              
                matelist.append(linealist)
                matenames.append(matename)
                linea = secfile.readline()


            elif linea[0:9] == '*LAMINATE':
                linea = linea.replace(" ", "")
                lamiid = lamiid + 1
                laminame = re.search('(?<=NAME=)\w+', linea).group()
                linea = secfile.readline()
                lamiangle = map(int, linea.split(','))
                if len(lamiangle) > maxlam:
                    maxlam = len(lamiangle)
                linea = secfile.readline()
                lamithick = map(float, linea.split(','))
                lamimate = secfile.readline()
                lamimate = lamimate.split(',')
                lamimate = [i.strip() for i in lamimate]
                lamidata = [laminame, lamiangle, lamithick, lamimate]
                lamilist.append(lamidata)
                laminames.append(laminame)
                linea = secfile.readline()

            else:
                print '- ERROR: Unidentified Keyword', linea,
                print 'Sorry, exiting...'
                break
        # -------------------------------------------------------------------------------------
        #     INITIALIZATION OF DICTIONARIES, ARRAYS AND LAMINATE OBJECTS
        # -------------------------------------------------------------------------------------
        # Form laminate objects
        lamiobj = []
        for i in np.arange(len(lamilist)):
            lamiobj.append( SectLami(lamilist[i]) )        
        self.lami = dict(zip(laminames,lamiobj))
        # Initializacion of material array
        mateobj = []
        for i in np.arange(len(matelist)):
            mateobj.append( SectMate(matelist[i]) )    
        self.mate = dict(zip(matenames,mateobj))
        # Initializacion of element array
        elsize = len(elemlist)
        eltype = [('el','intp',1),  ('n1','intp',1),  ('n2','intp',1),  ('type','a12',1), 
                  ('eset','a12',1), ('lset','a12',1), ('lath','float64',1),
                  ('inn1','float64',3), ('inn2','float64',3), 
                  ('mid1','float64',3), ('mid2','float64',3), 
                  ('out1','float64',3), ('out2','float64',3),
                  ('norm','float64',3), ('tang','float64',3), 
                  ('vec1','float64',3), ('vec2','float64',3)]
        self.elem = np.zeros(elsize, dtype=eltype)
        for i in np.arange(elsize):
            self.elem['el'][i] =  elemlist[i][0]
            self.elem['n1'][i] =  elemlist[i][1]
            self.elem['n2'][i] =  elemlist[i][2]
            self.elem['eset'][i] =  elemlist[i][3]  # Element Set Name
            self.elem['lset'][i] =  elemlist[i][4]  # Laminate Set Name
            self.elem['type'][i] =  elemlist[i][5]  # Laminate Set Name
        # Initializacion of node array
        ndsize = len(nodelist)
        ndtype = [('nd','intp',1),('coor','float64',3)]
        self.node = np.zeros(ndsize, dtype=ndtype)
        for i in np.arange(ndsize):
            self.node['nd'][i] = nodelist[i][0]
            self.node['coor'][i] = nodelist[i][1:4]
        print '-----> Input file was read in:', time.time() - start_time, 'seconds'

    def pregen(self):
        pregentime = time.time()
        ELEM = self.elem
        NODE = self.node
        # Reorder node and element tables
        old_nd_order = NODE['nd']
#        old_el_order = NODE['el']
        new_nd_order = np.arange(1, len(NODE) + 1)
        new_el_order = np.arange(1, len(ELEM) + 1) 
        NODE['nd'] = new_nd_order
        ELEM['el'] = new_el_order
        # Update Element table for node renumbering
        lookup_table = dict( zip( old_nd_order, new_nd_order ) ) # create your translation dict
        vect_lookup = np.vectorize( lookup_table.get ) # create a function to do the translation
        ELEM['n1'] = vect_lookup( ELEM['n1'] ) # Reassign the elements you want to change
        ELEM['n2'] = vect_lookup( ELEM['n2'] ) # Reassign the elements you want to change               
        # Correct Normals and Fill in element tangential a normal directions  
        numweb = len( np.extract(ELEM['type'] == 'WEB', ELEM['type'] ) )
        if numweb < 4:
            print '-----> ERROR: Minimum number of web elements per web is 2.'
        invcount = 0
        lyrcount = 0
        lstcount = []
        for el in ELEM:
            node1 = el['n1']
            node2 = el['n2']
            coor1 = NODE['coor'][node1 - 1]
            coor2 = NODE['coor'][node2 - 1]
            elnorm = ma.normal(coor1, coor2)           
            if np.dot( (0.5 * coor1 + 0.5 * coor2), elnorm.T) < 0.0:
#            if ma.dot( (0.5 * coor1 + 0.5 * coor2), elnorm) < 0.0:
                invcount += 1
                lstcount.append(el['el'])
                elnorm =  -elnorm
                el['n1'] = node2
                el['n2'] = node1               
                ctemp = coor1
                coor1 = coor2
                coor2 = ctemp
            el['norm'][:] = elnorm
            el['tang'][:] = ma.tangen(coor1, coor2)
            el['lath']= self.lami[el['lset']].lath
#            lyrcount += self.lami[el['lset']].lyrs['nlyr']
#        
#        # Initializacion of layers array        
#        self.lyrs = np.array
        print '-----> INFORMATION: Number of inverted normals corrected is:', invcount
        print '     Elements with inverted normal are:', lstcount
        print '-----> Pre was done in:', time.time() - pregentime, 'seconds'

    def amplif(self,ampy,ampz):
        """ Implifies the cross section """
        # Fastest version
        self.node['coor'] *=  np.array([0.0, ampy, ampz])
        # Slow version        
        #for i in np.arange(len(self.node)):
        #    coor_alias = self.node['coor'][i]
        #    coor_alias[1] *= ampy
        #    coor_alias[2] *= ampz

    def offset(self,offy,offz):
        """ Implifies the cross section """
        # Fastest version
        self.node['coor'] +=  np.array([0.0, offy, offz])
        # Fast version        
        #self.node['coor'][:,1] +=  offy
        #self.node['coor'][:,2] +=  offz
        # Slow version
        #for i in np.arange(len(self.node)):
        #    coor_alias = self.node['coor'][i]
        #    coor_alias[1] += offy
        #    coor_alias[2] += offz

    def writex(self, name):
        """ Write Output File """
        import csv
        nodedata = np.zeros( (len(self.node), 4) )
        nodedata[:,0] = self.node['nd']
        nodedata[:,1:4] = self.node['coor']
        elemdata = np.zeros( (len(self.elem), 3) )
        elemdata[:,0] = self.elem['el']
        elemdata[:,1] = self.elem['n1']
        elemdata[:,2] = self.elem['n2']

        f = open(name, "wb")
        # Node Output
        f.write('*NODE\r\n')
        for line in nodedata:
            csv.writer(f).writerow(line)
        # Element Output
        f.write('*ELEMENT,TYPE=T3D2\r\n')
        for line in elemdata:
            csv.writer(f).writerow(line)
        
        f.write('*END')
        print '-----> Mesh written to:', name

    def vectex(self):         
        """Calculationo of nodal normal vectors
        NOTE: Only el['out'] are neccessary, the inn and mid are redundant with layer"""
        ELEM = self.elem
        NODE = self.node
        iv = ma.unitx()
        jv = ma.unity()
        kv = ma.unitz()
        lyrcount = 0
        for el in ELEM:    
            n1 = el['n1']
            n2 = el['n2']
            coor1 = NODE['coor'][el['n1'] - 1]
            coor2 = NODE['coor'][el['n2'] - 1]   
            # WEB ELEMENTS FORMULATION            
            if el['type'] == 'WEB':
                el['vec1'][:] = el['norm']
                el['vec2'][:] = el['norm'] 
                factor = 1.5
                outproy1 = el['lath'] * el['vec1']
                outproy2 =  el['lath'] * el['vec2']
                halfoutproy1 = 0.5 * outproy1
                halfoutproy2 = 0.5 * outproy2               
                tangproy = factor * el['lath'] * el['tang']
                backidx = np.where(ELEM['n2'] == n1)[0]
                backel = ELEM['el'][backidx]
            #                    backcon = [j for j in backcon if j != elnum]
                forwidx = np.where(ELEM['n1'] == n2)[0]
                forwel = ELEM['el'][forwidx]
            #                    forwcon = [j for j in forwcon if j != elnum]
                if self.elem['type'][backidx] != 'WEB':
                    el['mid1'][:] =  coor1 + tangproy # Take on thickness for web joints
                    el['mid2'][:] =  coor2 
                    el['inn1'][:] =  el['mid1'] - halfoutproy1
                    el['inn2'][:] =  el['mid2'] - halfoutproy2
                    el['out1'][:] =  el['mid1'] + halfoutproy1
                    el['out2'][:] =  el['mid2'] + halfoutproy2
                if self.elem['type'][forwidx] != 'WEB' :
                    el['mid1'][:] =  coor1
                    el['mid2'][:] =  coor2 - tangproy
                    el['inn1'][:] =  el['mid1'] - halfoutproy1
                    el['inn2'][:] =  el['mid2'] - halfoutproy2
                    el['out1'][:] =  el['mid1'] + halfoutproy1
                    el['out2'][:] =  el['mid2'] + halfoutproy2 
            # RECTANGULAR ELEMENTS FORMULATION
            if el['type'] == 'R3D2':
                el['vec1'][:] = el['norm']
                el['vec2'][:] = el['norm']            
                # Proyections in the outer directions
                outproy1 = el['lath'] * el['vec1']
                outproy2 = el['lath'] * el['vec2']
                halfoutproy1 = 0.5 * outproy1
                halfoutproy2 = 0.5 * outproy2
                el['inn1'][:] =  coor1 - outproy1
                el['inn2'][:] =  coor2 - outproy2
                el['mid1'][:] =  coor1 - halfoutproy1
                el['mid2'][:] =  coor2 - halfoutproy2
                el['out1'][:] =  coor1
                el['out2'][:] =  coor2                           
            # FORMULATION FOR TRAPEZOIDAL ELEMENTS
            if el['type'] == 'T3D2':
                backidx = np.where(ELEM['n2'] == n1)[0]
                backel = ELEM['el'][backidx]
                forwidx = np.where(ELEM['n1'] == n2)[0]
                forwel = ELEM['el'][forwidx]
                thckel = el['lath']                   
                # Backward and Forward element nodes coordinates
                coorb1 = NODE['coor'][ELEM['n1'][backidx] - 1][0]
                coorb2 = NODE['coor'][ELEM['n2'][backidx] - 1][0]
                coorf1 = NODE['coor'][ELEM['n1'][forwidx] - 1][0]
                coorf2 = NODE['coor'][ELEM['n2'][forwidx] - 1][0]                         
                # COVARIANT COMPONENT
                covt = np.array( [0.0, thckel, thckel] )
                # COLINEAR VECTORS (BACKWARD AND FORWARD)
                gsubeb = coorb2 - coorb1
                gsubb  = coor1  - coor2
                gsubf  = coor2  - coor1 
                gsubef = coorf1 - coorf2 
                # NORMALIZED PERPENDICULAR VECTORS
                gsupb  = np.array( [ 0.0, gsubeb[2], -gsubeb[1] ] )   
                gsupb  = gsupb / ma.nrm(gsupb)                   
                gsupeb = -1.0 * np.array( [ 0.0, gsubb[2],  -gsubb[1]  ] )    
                gsupeb = gsupeb/ma.nrm(gsupeb)        
                gsupf  = -1.0 * np.array( [ 0.0, gsubef[2], -gsubef[1] ] )  
                gsupf  = gsupf/ma.nrm(gsupf)
                gsupef = np.array( [ 0.0,  gsubf[2], -gsubf[1] ] )   
                gsupef = gsupef/ma.nrm(gsupef)     
                axial  = ma.unitx()    
                # RECIPROCAL BASIS
                vb = np.dot( gsupeb, ma.crs(axial, gsupb) )  # Volume or paralellepiped
                vf = np.dot( gsupf,  ma.crs(axial, gsupef) )  # Volume or paralellepiped          
                # COLINEAR NODES BYPASS (ZERO PARALLELEPIDED VOLUME)
                vbtol, vftol = 1E-3, 1E-3 
                vbtolcount, vftolcount = 0, 0                 
                if vb <= vbtol:
#                    print '-----> WARNING!. Negative parallelepiped volume for backward element:', el['el']
                    vbtolcount += 1
                    el['vec1'][:] = ma.versor( 0.5 * (gsupb + gsupeb) )
                    outproy1 = el['lath'] * el['vec1']
                    el['out1'][:] =  coor1                        
                    el['mid1'][:] =  coor1 - 0.5 * outproy1
                    el['inn1'][:] =  el['out1'] - outproy1
                else:
                    # Reciprocal Bases    
                    eb1 = ma.crs(gsupb, gsupeb) / vb  
                    eb2 = ma.crs(gsupeb, axial) / vb  
                    eb3 = ma.crs(axial,  gsupb) / vb                         
                    # TRANSFORMATION
                    Mb = np.array ( [[ np.dot(iv, eb1),  np.dot(jv, eb1),  np.dot(kv, eb1)],
                                     [ np.dot(iv, eb2),  np.dot(jv, eb2),  np.dot(kv, eb2)],
                                     [ np.dot(iv, eb3),  np.dot(jv, eb3),  np.dot(kv, eb3)]])  # Transformation matrix, Original-Cartesian
                    tvb = np.dot( Mb.T, covt)  # Cartesian components of vector 
                    # NEW POSITION OF ELEMENT NODES
                    el['vec1'][:] = tvb / thckel
                    el['out1'][:] = coor1                        
                    el['mid1'][:] = coor1 - 0.5 * tvb
                    el['inn1'][:] = coor1 - tvb
                if vf <= vftol:
#                    print '-----> WARNING!. Negative parallelepiped volume for forward element :', el['el']
                    vftolcount += 1          
                    el['vec2'][:] = ma.versor( 0.5 * (gsupf + gsupef) )
                    outproy2 = el['lath'] * el['vec2']
                    el['out2'][:] =  coor2
                    el['mid2'][:] =  coor2 - 0.5 * outproy2
                    el['inn2'][:] =  el['out2'] - outproy2
                else:
                    # Reciprocal Bases    
                    ef1 = ma.crs(gsupef, gsupf) / vf
                    ef2 = ma.crs(gsupf, axial) / vf
                    ef3 = ma.crs(axial,  gsupef) / vf
                    # TRANSFORMATION
                    Mb = np.array ( [[ np.dot(iv, ef1),  np.dot(jv, ef1),  np.dot(kv, ef1) ],
                                     [ np.dot(iv, ef2),  np.dot(jv, ef2),  np.dot(kv, ef2) ],
                                     [ np.dot(iv, ef3),  np.dot(jv, ef3),  np.dot(kv, ef3) ]])  # Transformation matrix, Original-Cartesian
                    tvf = np.dot( Mb.T, covt)  # Cartesian components of vector        
                    # NEW POSITION OF ELEMENT NODES
                    el['vec2'][:] = tvf / thckel # NOTE: This is not a versor, this is a thickness unit vector
                    el['out2'][:] = coor2 
                    el['mid2'][:] = coor2 - 0.5 * tvf
                    el['inn2'][:] = coor2 - tvf

    def layers(self):
        lyrcount = 0
        jl = 0    
        ELEM = self.elem
        for el in ELEM:
            lyrcount += self.lami[el['lset']].nlyr
        lyrtype = [ ('nlyr','intp',1), ('nele','intp',1), 
                    ('mate','a12',1),  ('angl','float64',1),
                    ('inn1','float64',3), ('inn2','float64',3),
                    ('mid1','float64',3), ('mid2','float64',3),
                    ('mad1','float64',3), ('mad2','float64',3),
                    ('out1','float64',3), ('out2','float64',3),
                    ('long','float64',1), ('thck','float64',1),
                    ('ntop','float64',1), ('nmid','float64',1), ('nbot','float64',1),
                    ('norm','float64',3), ('tang','float64',3),
                    ('rn_i','float64',1), ('rn_s','float64',1),
                    ('rs_1','float64',1), ('rs_2','float64',1)]
        self.lyrs = np.zeros(lyrcount, dtype=lyrtype)        
        LYRS = self.lyrs
        for elid, el in enumerate(ELEM):
            lami = self.lami[el['lset']]
            for ly in lami.lyrs:
                LYRS['mate'][jl] = ly['mate']
                LYRS['thck'][jl] = ly['thck']
                LYRS['angl'][jl] = ly['angl']
                
                LYRS['out1'][jl] = el['out1'] + ly['ntop'] * el['vec1']
                LYRS['mid1'][jl] = el['out1'] + (ly['ntop'] - 0.5 * ly['thck'])  * el['vec1']
                LYRS['mad1'][jl] = el['out1'] + (ly['ntop'] - 0.45 * ly['thck']) * el['vec1']                
                LYRS['inn1'][jl] = el['out1'] + ly['nbot'] * el['vec1']
                 
                LYRS['out2'][jl] = el['out2'] + ly['ntop'] * el['vec2']
                LYRS['mid2'][jl] = el['out2'] + (ly['ntop'] - 0.5 * ly['thck'])  * el['vec2']
                LYRS['mad2'][jl] = el['out2'] + (ly['ntop'] - 0.45 * ly['thck']) * el['vec2']                  
                LYRS['inn2'][jl] = el['out2'] + ly['nbot'] * el['vec2']
                
                LYRS['long'][jl] = ma.nrm( LYRS['mid2'][jl] - LYRS['mid1'][jl])
                
                LYRS['ntop'][jl] = ly['ntop']
                LYRS['nmid'][jl] = ly['nmid']
                LYRS['nbot'][jl] = ly['nbot']
                
                LYRS['norm'][jl] = el['norm']
                LYRS['tang'][jl] = el['tang']                
                
                LYRS['rn_i'][jl] = ma.dot( el['norm'], LYRS['inn1'][jl] ) 
                LYRS['rn_s'][jl] = ma.dot( el['norm'], LYRS['out1'][jl] ) 
                LYRS['rs_1'][jl] = ma.dot( el['tang'], LYRS['mid1'][jl] ) 
                LYRS['rs_2'][jl] = ma.dot( el['tang'], LYRS['mid2'][jl] ) 
                
                jl += 1
                LYRS['nlyr'][jl - 1] = jl
                LYRS['nele'][jl - 1] = elid
            lyrcount += self.lami[el['lset']].nlyr
        
    def center(self): 
        """ Centers calculation """
        ELEM = self.elem
        LYRS = self.lyrs
        MATE = self.mate
        centype = [('area', 'intp', 1), ('mass', 'float64', 1), ('caxs', 'intp', 3), ('cgxs', 'float64', 3)]
        self.cent = np.zeros(1, dtype=centype)        
        mass = 0.0
        area = 0.0
        cgxs = np.zeros(3)
        caxs = np.zeros(3)

        for ly in LYRS:
            arealy = ly['thck'] * ly['long']
            massly = arealy * MATE[ly['mate']].dens
            coorly = 0.5 * (ly['mid1'] + ly['mid2'])
            area += arealy
            mass += massly
            caxs += coorly * arealy
            cgxs += coorly * massly
        self.cent['area'] = area
        self.cent['mass'] = mass
        self.cent['caxs'] = caxs / area
        self.cent['cgxs'] = cgxs / mass
        print '-----> INFORMATION: Total mass of the cross section is:', mass
        print '-----> INFORMATION: Total area of the cross section is:', area
        print '-----> INFORMATION: Center of mass of the cross section is:', self.cent['cgxs'] 
        print '-----> INFORMATION: Center of area of the cross section is:', self.cent['caxs'] 

    def stiff(self):
        DSEC = np.zeros((6, 6))
        DLYR = np.zeros((6, 6))
        JSEC = np.zeros((3, 3))
        MSEC = np.zeros((3, 3))
        LYRS = self.lyrs
        MATE = self.mate
        shrc = 0.0
        self.DLYR = np.zeros((6, 6, len(LYRS)))
        for ly in LYRS:
            E1  = MATE[ly['mate']].elas[0]
            E3  = MATE[ly['mate']].elas[1]
            G12 = MATE[ly['mate']].elas[2]
            G13 = G12
            G23 = MATE[ly['mate']].elas[3]
            v13 = MATE[ly['mate']].elas[4]
            rho = MATE[ly['mate']].dens
            mm = np.cos( ly['angl'] * np.pi / 180.0 )
            nn = np.sin( ly['angl'] * np.pi / 180.0 )

            # COMPLIANCE PARA ESTADO PLANO DE TENSIONES EN LA LAMINA
            SL = np.array( [[   1/E1,      -v13/E1,      0    ],
                            [  -v13/E1,       1/E3,      0    ],
                            [    0,            0,     1/G13   ]]) 
    
            # Transforamcion de tensiones
            TS = np.array( [[  mm**2,   nn**2,     2*mm*nn    ],
                            [  nn**2,   mm**2,    -2*mm*nn    ],
                            [-mm*nn,   mm*nn,    mm**2-nn**2  ]])
    
            TE = np.array( [[  mm**2,      nn**2,      mm*nn    ],
                            [  nn**2,      mm**2,     -mm*nn    ],
                            [-2*mm*nn,   2*mm*nn,   mm**2-nn**2 ]])   

            SS = np.dot( np.dot(TE, SL), np.linalg.inv(TS) )  # Cuidado Mogotron! Use inverse and not transpose since the transformation IS NOT ORTHOGONAL

            # Estado uniaxial de tensiones, o sea tension circunferencial cero    
            SSU = np.array([[ SS[0,0], SS[0,2] ], 
                            [ SS[2,0], SS[2,2] ]])       

            CL = np.linalg.inv( SSU )        

            C11 = CL[0,0]  
            C12 = CL[0,1]  
            C21 = CL[1,0] 
            C22 = CL[1,1]  

            # Transformation Tensor
            Q = np.array( [ ma.unitx(), ly['norm'], ly['tang'] ] ).T

            n0, n1 = ly['nbot'], ly['ntop'] # n axis coordinates      
            rni, rns = ly['rn_i'],  ly['rn_s']   # Inner and Outer Projections in the normal direction
            rs1, rs2 = ly['rs_1'],  ly['rs_2']

            # No warping constants
            c1 = (rns-rni)*(rs2-rs1);
            c2n = 0.5*(rns**2-rni**2)*(rs2-rs1);
            c2s = 0.5*(rns-rni)*(rs2**2-rs1**2);
            c3n = (1.0/3.0)*(rns**3-rni**3)*(rs2-rs1);
            c3s = (1.0/3.0)*(rns-rni)*(rs2**3-rs1**3);
            c4 = 0.25*(rns**2-rni**2)*(rs2**2-rs1**2);

            # Constants with warping
            dws, cwarp, cwcoup = 0, 0, 0  # Warping and Warping Couplings FLAGS

            crnw = dws * (rns**2-rni**2)*(rs2-rs1) + dws**2 * (rns-rni)*(rs2-rs1);
            cn2w = cwarp * 0.5 *dws*(rns**2-rni**2)*(rs2-rs1);
            cs2w = cwarp * 0.5 *dws*(rns-rni)*(rs2**2-rs1**2);
            c1w =  cwarp * c1* dws;

            D11 = np.array( [[   C11*c1,       0.0,      C12*c1 ],
                             [   0.0,          0.0,      0.0    ],           
                             [   C12*c1,       0.0,      C22*c1 ]])

            D22 = np.array( [[  C22*(c3n+cwarp*crnw),   C12*(c3n+cwarp*cn2w),   C12*(c4+cwarp*cs2w)], 
                             [  C12*(c3n+cwarp*cn2w),   C11*c3n,                C11*c4             ],     
                             [  C12*(c4+cwarp*cs2w),    C11*c4,                 C11*c3s            ]])  

            D12 = np.array( [[  C12*(c2n+cwcoup*c1w),   C11*c2n,    C11*c2s], 
                             [  0.0,                    0.0,        0.0    ],           
                             [  C22*(c2n+cwcoup*c1w),   C12*c2n,    C12*c2s]])

            
            
            DLYR[0:3,0:3] = ma.do3(Q, D11, Q.T)
            DLYR[3:6,3:6] = ma.do3(Q, D22, Q.T)
            DLYR[0:3,3:6] = ma.do3(Q, D12, Q.T)
            DLYR[3:6,0:3] =  DLYR[0:3,3:6].T 
            self.DLYR[0:6,0:6, ly['nlyr']-1] = DLYR / (ly['thck'] * ly['long']) # Stiffness Density
            
            DSEC += DLYR

            # INERTIA MATRIX
            jxs = np.array ([[c3n+c3s,    0.0,   0.0  ],
                             [0.0,        c3n,   -c4  ],
                             [0.0,        -c4,   c3s  ]])

            JSEC += ma.do3(Q.T, rho * jxs, Q)
            MSEC += ma.dia( ly['thck'] * ly['long'] * rho)

            shrc += C11 * c1 * 0.5 * ( ly['mid1'] + ly['mid2'] ) 
            
    
        self.ISEC = np.zeros((6, 6))
        self.ISEC[0:3, 0:3] = MSEC
        self.ISEC[3:6, 3:6] = JSEC        
        self.DSEC = DSEC
        self.shrc = shrc / DSEC[0, 0]

    def plotxs(self,ptype, *pdof):
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        import matplotlib as mpl
        ax = plt.subplot()
        ELEM = self.elem
        vert = np.zeros((len(ELEM), 4, 2))
        numpoly = len(ELEM)
        # Plot Midsurface
        if ptype == 'mids':
            ELEM = self.elem
            vert = np.zeros((len(ELEM), 4, 2))
            numpoly = len(ELEM)
            for i, el in enumerate(ELEM):
                vert[i,0,:] = el['mid1'][1:3]
                vert[i,1,:] = el['mid1'][1:3] + 0.002 * el['norm'][1:3]  
                vert[i,2,:] = el['mid2'][1:3] + 0.002 * el['norm'][1:3] 
                vert[i,3,:] = el['mid2'][1:3]
            color = np.random.random(numpoly) * 500 
        
        # Plot geometry
        if ptype == 'geom':
            ELEM = self.elem
            vert = np.zeros((len(ELEM), 4, 2))
            numpoly = len(ELEM)
            for i, el in enumerate(ELEM):          
                vert[i,0,:] = el['out1'][1:3]
                vert[i,1,:] = el['inn1'][1:3]
                vert[i,2,:] = el['inn2'][1:3]
                vert[i,3,:] = el['out2'][1:3]
            color = np.random.random(numpoly) * 500 
        
        # Plot Layer by color    
        if ptype == 'lyrs' or ptype == 'mlyr' or ptype == 'stif':           
            matcolor = np.arange(len(self.mate)) + 1
            matecolor = dict(zip(self.mate.keys(), matcolor))
            LYRS = self.lyrs
            vert = np.zeros((len(LYRS), 4, 2))
            color = np.zeros(len(LYRS))
            numpoly = len(LYRS) 
            
            if ptype == 'lyrs' or ptype == 'stif':
                for i, el in enumerate(LYRS):          
                    vert[i,0,:] = el['out1'][1:3]
                    vert[i,1,:] = el['inn1'][1:3]
                    vert[i,2,:] = el['inn2'][1:3]
                    vert[i,3,:] = el['out2'][1:3]         
                    if ptype == 'lyrs':
                        color[i] = matecolor[el['mate']] #Color by mat 
                        matecolor[el['mate']] 
                        ax.set_title('Geometry Material Mapping Contour')
                    if ptype == 'stif':
                        color[i] = self.DLYR[pdof[0]-1, pdof[0]-1, i]  #Color by stiffness
                        ax.set_title('Layer Stiffness Density Contour in Direction: ' + str(pdof[0]))
            elif ptype == 'mlyr':
                for i, el in enumerate(LYRS):
                    vert[i,0,:] = el['mid1'][1:3]
                    vert[i,1,:] = el['mad1'][1:3]   
                    vert[i,2,:] = el['mad2'][1:3]  
                    vert[i,3,:] = el['mid2'][1:3]
                    color[i] = 0
        # Execute Plot     
#        ax = plt.subplot()
        coll = PolyCollection(vert, array=color, cmap=mpl.cm.jet, edgecolors='black')
        coll.set_array(np.array(color))
        ax.add_collection(coll)
        ax.autoscale_view()
        ax.axis('equal')        
        plt.colorbar(coll)
            
class SectLami(object):
    def __init__(self,datalist):
        lyrnum = len(datalist[1])
        self.name = datalist[0]
        self.lath = sum(datalist[2])
        self.nlyr = lyrnum
        self.mate = datalist[3]
        # Organize lamination angles, thicknesses and coordinates in numpy array
        lyrtype = [('nlyr','intp',1), ('mate','a12',1), ('angl','float64',1),  
                   ('thck','float64',1), ('nbot','float64',1), ('nmid','float64',1), ('ntop','float64',1)]
        self.lyrs = np.zeros(lyrnum, dtype=lyrtype)
        self.lyrs['nlyr'] = np.arange(1, lyrnum+1)
        self.lyrs['mate'] = datalist[3]
        self.lyrs['angl'] = np.array( datalist[1] )
        self.lyrs['thck'] = np.array( datalist[2] )
        self.lyrs['nbot'] = -1.0 * np.cumsum( datalist[2] )
        self.lyrs['nmid'] = self.lyrs['nbot'] + 0.5 * self.lyrs['thck']
        self.lyrs['ntop'] = self.lyrs['nbot'] + self.lyrs['thck']
##    def __repr__(self,datalist):
#        return str(('Object ', self.name))
        
class SectMate(object):
    def __init__(self,datalist):
        self.name = datalist[0]
        self.type = datalist[1]
        if self.type == 'ISOTROPIC':
            E1 = datalist[2][0]
            v13 = datalist[2][1]
            G13 = E1 / (2 + 2 * v13)
            self.elas = [E1, E1, G13, G13, v13 ]
            self.dens = datalist[2][2]
        if self.type == 'TRISOTROPIC':
            self.elas = datalist[2][0:5]
            self.dens = datalist[2][5]
#    def __repr__(self,datalist):
#        return str(('Object:', self.name))