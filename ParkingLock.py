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
        c.MPLPlot(style='--ob')
        

    
class Gear:
    
    def __init__(self,Ox,Oy,Z,DI,DF,L8,L9,theta1,theta2,theta3,alphaJ,alphaR):
        self.Ox = Ox
        self.Oy = Oy
        self.Z = Z
        self.DI = DI
        self.DF = DF
        self.L8 = L8
        self.L9 = L9 #Distance JK
        self.theta1 = theta1
        self.theta2 = theta2 
        self.theta3 = theta3 #(theta2+theta3)=(2*npy.pi/Z) ?
        self.alphaJ = alphaJ
        self.aplhaR = alphaR
        
        
        
       
    def Geometry(self):
        theta = 2*npy.pi/Z #theta=theta2 +theta3
        
        R = O.Translation(vm.point2D((self.L8*npy.cos(alphaR),-self.L8*npy.sin(alphaR))))
        
        K = R.Translation(vm.point2D((0,DF/2)))
        H = K.Rotation(R,theta2)
        M = K.Rotation(R,theta2+theta3)
        a1 = vm.Arc2D(K,H,M)
        
        J1 = R.Translation(vm.Point2D((0,DI/2)))
        J=J1.Rotation(R,(theta2-theta1)/2)
        #est ce un probleme de ne pas controler directement la pente du créneau ? (on ne se sert pas de alphaJ)
        #sinon on pourrait remplacer les deux lignes précédentes par:J=K.Translation(vm.point2D((-self.L9*npy.sin(alphaJ),-self.L9*npy.cos(alphaJ))))
        #Dans ce cas il faudrait peut être paramétrer H de la même façon
        I=J.Rotation(R,theta1)
        
        l1=primitives2D.RoundedLineSegments2D([K,J,I,H],{'''je ne sais pas trop quoi mettre ici'''})
        
        L=[a1,l1]
        for i in range(Z-1):
            thetar=(i+1)*theta
            L.append(a1.Rotation(R,thetar,True))
            L.append(l1.Rotation(R,thetar,True))

        c1=vm.Contour2D(L)
        c1.MPLPlot()

        
P1 = Pawl(Ox = 0, Oy = 0, alpha0 = npy.pi/4., LA = 1,
          alphaA = npy.pi/4., L1 = 0.1, alphaB = npy.pi/4., L2 = 0.1,
          alphaP = npy.pi/6., L5 = 0.1, alphaD = npy.pi/6., LE = 0.8, 
          alpha1 = npy.pi/4. - npy.pi/8., L3 = 0.05, L4 = 0.05, alphaE = npy.pi/8., L6 = 0.1, 
          alphaC = npy.pi/6., L7 = 0.1, DF = 1, DI = 0.9)
P1.Geometry()

#G1 = Gear(teeth_number = 12, pitch_diameter = 0.07, point_o = [0, 0], 
#          external_diameter = 0.08, internal_diameter = 0.06)

#PL1 = ParlinkLock(pawl = P1, point_a = [0, 1], alpha = 0.2, gear = G1, point_o = [0, 0], betha = 0)
