# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""
import volmdlr.primitives2D as primitives2D
from scipy.optimize import minimize, fsolve
import numpy as npy
import matplotlib.pyplot as plt

import volmdlr as vm
import volmdlr.primitives2D as primitives2D

class Pawl:
    """
    point_a : centre de rotation du doigt
    alpha0 : angle entre x et le plan du plateau
    """
    def __init__(self, Ox, Oy, alpha0, LA, alphaA, L1, alphaB, L2,
                 alphaP, L5, alphaD, LE, alpha1, L3, L4, alphaE, L6, 
                 alphaC, L7, DF, DI):
        self.Ox = Ox
        self.Oy = Oy
        self.alpha0 = alpha0
        self.LA = LA
        self.alphaA = alphaA
        self.L1 = L1
        self.alphaB = alphaB
        self.L2 = L2
        self.alphaP = alphaP
        self.L5 = L5
        self.alphaD = alphaD
        self.LE = LE
        self.alpha1 = alpha1
        self.L3 = L3
        self.L4 = L4
        self.alphaE = alphaE
        self.L6  = L6
        self.alphaC = alphaC
        self.L7 = L7
        self.DF = DF
        self.DI = DI
        
    def Geometry(self):
        
        O  = vm.Point2D((self.Ox, self.Oy))
        A  = O.Translation(vm.Point2D((self.LA*npy.cos(self.alpha0), self.LA*npy.sin(self.alpha0))))
        alpha = self.alpha0 - self.alphaA
        B  = A.Translation(vm.Point2D((self.L1*npy.cos(alpha), self.L1*npy.sin(alpha))))
        alpha = alpha - self.alphaB
        P  = B.Translation(vm.Point2D((self.L2*npy.cos(alpha), self.L2*npy.sin(alpha))))
        alpha = alpha + self.alphaP
        M  = P.Translation(vm.Point2D((self.L5*npy.cos(alpha), self.L5*npy.sin(alpha))))
        
        alpha = self.alpha0 - self.alphaE
        E  = O.Translation(vm.Point2D((self.LE*npy.cos(alpha), self.LE*npy.sin(alpha))))
        alphap = alpha - self.alpha1
        C  = E.Translation(vm.Point2D((self.L4*npy.cos(alphap), self.L4*npy.sin(alphap))))
        D  = E.Translation(vm.Point2D((-self.L3*npy.cos(alphap), -self.L3*npy.sin(alphap))))
        alpha = alphap + npy.pi/2. - self.alphaC
        F  = C.Translation(vm.Point2D((self.L7*npy.cos(alpha), self.L7*npy.sin(alpha))))
        alpha = alphap + npy.pi/2. + self.alphaD
        G  = D.Translation(vm.Point2D((self.L6*npy.cos(alpha), self.L6*npy.sin(alpha))))
        
        AB=vm.LineSegment2D(A, B)
        BP=vm.LineSegment2D(B, P)
        PM=vm.LineSegment2D(P, M)
        DE=vm.LineSegment2D(D, E)
        EC=vm.LineSegment2D(E, C)
        DG=vm.LineSegment2D(D, G)
        CF=vm.LineSegment2D(C, F)
        
        alpha = alphap - npy.pi/2.
        R = E.Translation(vm.Point2D((self.DI/2.*npy.cos(alpha), self.DI/2.*npy.sin(alpha))))
        # recherche du point F1 et F2
        def fonct(alpha, *data):
            line, d = data
            x  = R.Translation(vm.Point2D((d/2.*npy.cos(alpha + npy.pi/2.), d/2.*npy.sin(alpha + npy.pi/2.))))
            p = line.PointProjection(x)
            return p.PointDistance(x)
        sol1 = fsolve(fonct, 0, args=(DG, self.DF))
        F1  = R.Translation(vm.Point2D((self.DF/2.*npy.cos(npy.pi/2. + sol1), self.DF/2.*npy.sin(npy.pi/2. + sol1))))
        sol2 = fsolve(fonct, 0, args=(CF, self.DF))
        F2  = R.Translation(vm.Point2D((self.DF/2.*npy.cos(npy.pi/2. + sol2), self.DF/2.*npy.sin(npy.pi/2. + sol2))))
        
        Cr = vm.Circle2D(R, self.DF/2.)
        
        c=vm.CompositePrimitive2D([O, AB, BP, PM, DE, EC, DG, CF, R, Cr, F1, F2])
        return c
    
    def Plot(self):
        c = self.Geometry()
        c.MPLPlot(style='--ob')
        

    
class Gear:
    
    def __init__(self, Ox, Oy, Z, DI, DF, DE, 
                 alpha1, alpha2, alpha3):
        self.Ox = Ox
        self.Oy = Oy
        self.Z = Z
        self.DI = DI
        self.DF = DF
        self.DE = DE
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.alpha3 = alpha3
        
       
    def Geometry(self):
        
        pas = 2*npy.pi/self.Z
        O = vm.Point2D((self.Ox, self.Oy))
        A = O.Translation(vm.Point2D((0, self.DE/2.)))
        Bp = O.Translation(vm.Point2D((0, self.DI/2.)))
        B = Bp.Rotation(O, -self.alpha1, copy = True)
        Bi = B.Rotation(O, -self.alpha2/2., copy = True)
        C = B.Rotation(O, -self.alpha2, copy = True)
        D = A.Rotation(O, -self.alpha1 - self.alpha2 - self.alpha3, copy = True)
        Di = D.Rotation(O, (-pas + self.alpha1 + self.alpha2 + self.alpha3)/2., copy = True)
        Ap = A.Rotation(O, -pas, copy = True)
        
        dir1 = (B[1] - A[1])/(B[0] - A[0])
        ord1 = B[1] - dir1*B[0]
        def fonct(x, *data):
            Ox, Oy, R, coeff_dir, ord_origine = data
            f = (ord_origine + coeff_dir*x - Oy)**2 + (x - Ox)**2 - R**2
            return f
        x1 = fsolve(fonct, 0, args=(self.Ox, self.Oy, self.DF/2., dir1, ord1))
        R1 = vm.Point2D((x1, ord1 + x1*dir1))
        
        dir2 = (D[1] - C[1])/(D[0] - C[0])
        ord2 = C[1] - dir2*C[0]
        x2 = fsolve(fonct, 1, args=(self.Ox, self.Oy, self.DF/2., dir2, ord2))
        R2 = vm.Point2D((x2, ord2 + x2*dir2))
        
        li = [B, Bi, C, D, Di, Ap]
        AB=vm.LineSegment2D(A, B)
        BC = vm.Arc2D(B, Bi, C)
        CD=vm.LineSegment2D(C, D)
        
        list_seg = []
        Am = A
        for z in range(self.Z):
            li_point = []
            for pt in li:
                li_point.append(pt.Rotation(O, -pas*z, copy=True))
            AB = vm.LineSegment2D(Am, li_point[0])
            BC = vm.Arc2D(li_point[0], li_point[1], li_point[2])
            CD=vm.LineSegment2D(li_point[2], li_point[3])
            DA = vm.Arc2D(li_point[3], li_point[4], li_point[5])
            list_seg.extend([AB, BC, CD, DA])
            Am = li_point[5]
                
        c=vm.CompositePrimitive2D(list_seg)
        return c
    
    def Plot(self):
        c = self.Geometry()
        c.MPLPlot(style='--ob')

class ParkingLock:
    def __init__(self, pawl, gear, pos_pawl, alpha_pawl, pos_gear, alpha_gear):
        
        self.pawl = pawl
        self.gear = gear
        self.pos_pawl = pos_pawl
        self.alpha_pawl = alpha_pawl
        self.pos_gear = pos_gear
        self.alpha_gear = alpha_gear
    
#    def Geometry(self):
        
        
        
P1 = Pawl(Ox = 0, Oy = 0, alpha0 = npy.pi/4., LA = 1,
          alphaA = npy.pi/4., L1 = 0.1, alphaB = npy.pi/4., L2 = 0.1,
          alphaP = npy.pi/6., L5 = 0.1, alphaD = npy.pi/6., LE = 0.8, 
          alpha1 = npy.pi/4. - npy.pi/8., L3 = 0.05, L4 = 0.05, alphaE = npy.pi/8., L6 = 0.1, 
          alphaC = npy.pi/6., L7 = 0.1, DF = 1, DI = 0.9)
P1.Plot()

G1 = Gear(Ox = 0, Oy = 0, Z = 10, DI = 0.1, DF = 0.11, DE = 0.12, 
          alpha1 = 0.1, alpha2 = 0.3, alpha3 = 0.1)
G1.Plot()

#G1 = Gear(teeth_number = 12, pitch_diameter = 0.07, point_o = [0, 0], 
#          external_diameter = 0.08, internal_diameter = 0.06)

#PL1 = ParlinkLock(pawl = P1, point_a = [0, 1], alpha = 0.2, gear = G1, point_o = [0, 0], betha = 0)
