from fractions import Fraction
def getMatrixDet(matrix):
    if len(matrix) == 1:
        return matrix[0][0]
    if len(matrix) == 2:
        return matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]
    detMatrix = 0
    for i in range(len(matrix)):
        minorMatrix = getMinor(matrix, 0, i)
        minorDet = getMatrixDet(minorMatrix)
        detMatrix += ((-1)**(i))*float(matrix[0][i])*(minorDet)
    return detMatrix

def getTranspose(matrix):
    matrixTranspose = []
    for col in range(len(matrix[0])):
        transposeRow = []
        for row in range(len(matrix)):
            transposeRow.append(matrix[row][col])
        matrixTranspose.append(transposeRow)
    return matrixTranspose        
		
def getMinor(matrix, i, j):
    minorMatrix = []
    for row in range(len(matrix)):
        minorRow = []
        for col in range(len(matrix)):
            if row == i or col==j:
                continue
            minorRow.append(matrix[row][col])
        if len(minorRow)!=0:
            minorMatrix.append(minorRow)
    return minorMatrix	

def getFirstLineOfMatrixInverse(detMatrix, matrix):
    firstRow = []
    minorsStack = []
    for i in range(len(matrix)):
        minorsStack.append(getMatrixDet(getMinor(matrix, i, 0)))
    for i in range(len(minorsStack)):
        minor = minorsStack[i]
        firstRow.append(((-1)**i)*minor*(float(1)/detMatrix))
    return firstRow		

def lcm(a, b):
    return abs(a*b) // gcd(a, b)

def gcd(a, b):
    if (b == 0):
        return a
    return gcd(b, a % b)

def diffMatrix(A, B):
    if len(A) != len(B) or len(A[0]) != len(B[0]):
        return []
    diffAB = []
    for irow in range(len(A)):
        diffABRow = []
        for jcol in range(len(A[0])):
            diffABRow.append(A[irow][jcol]-B[irow][jcol])
        diffAB.append(diffABRow)
    return diffAB
         
def getProbMatrix(matrix):
    probMatrix = []
    for row in matrix:
        probMatrixRow = []
        sumOfWeights = sum(row)
        if sumOfWeights != 0:
            for w in row:
                probMatrixRow.append(float(w)/sumOfWeights)
            probMatrix.append(probMatrixRow) 
        else:
            probMatrix.append(row)        
    return probMatrix

def ones(n):
    ones = []
    for i in range(n):
        onesRow = [0]*n
        onesRow[i] = 1
        ones.append(onesRow)
    return ones
    
def splitMatrixIntoQR(matrix):
    Q = []
    R = []
    absorbing = set()
    for i in range(len(matrix)):
        sumOfWeights = sum(matrix[i])
        if sumOfWeights == 0:
            absorbing.add(i)
    for i in range(len(matrix)):
        #print(i)
        if i in absorbing:
             continue
        sumOfWeights = float(sum(matrix[i]))
        QRow = []
        RRow = []
        for j in range(len(matrix[0])):
            if j in absorbing:
                RRow.append(float(matrix[i][j]/sumOfWeights))
            else:
                QRow.append(float(matrix[i][j]/sumOfWeights))
        Q.append(QRow)
        R.append(RRow)
    return Q, R

def solution(m):
    # *** Markov Chains Absorbing ***
    if len(m) < 2:
        return [1,1]
    Q, R = splitMatrixIntoQR(m)
    I = []
    I = ones(len(Q)) 
    #print(Q)
    #print(R)
    #print(I)
    Rt = getTranspose(R)
    N = diffMatrix(I, Q)
    detN = getMatrixDet(N)
    firstRow = getFirstLineOfMatrixInverse(detN, N)
    decProb = []
    fracRes = []
    sameDen = 1
    for j in range(len(Rt)):
        cofact = 0
        for l in range(len(firstRow)):
            cofact += firstRow[l]*Rt[j][l]
        decProb.append(Fraction(cofact).limit_denominator())
    
    for i in range(len(decProb)-1):
        curDen = lcm(decProb[i].denominator, decProb[i+1].denominator)
        if sameDen < curDen:
            sameDen = curDen
    for r in decProb:
        fracRes.append(int(r.numerator*(sameDen//r.denominator)))
    fracRes.append(int(sameDen))    
    return fracRes
