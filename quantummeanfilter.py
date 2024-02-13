import classiq
from classiq import Model,synthesize,show,execute,set_constraints
from classiq.builtin_functions import HGate,CXGate,IGate,Mcx,XGate,SwapGate,Adder,CCXGate,ComputationalBasisStatePreparation
from classiq.quantum_register import QReg,QUInt
from classiq.interface.generator.mcu import CTRL,TARGET
from classiq.model import Preferences,Constraints,OptimizationParameter
from classiq.execution import (
    ExecutionPreferences,
    ExecutionDetails,
    IBMBackendPreferences,
)
from PIL import Image
import numpy as np
import math

class ImageToArray:
     def getImage(self,path,n:int=2):
      img=Image.open(path)# image path
      new_size = (n, n)
      gray_img = img.convert("L")
      resized_image=gray_img.resize(new_size)
      _image=np.array(resized_image)
      return _image


class ImageToCircuit:
   __toBinary= lambda self,x,size:bin(x)[2:].zfill(size)
   __applyLog= lambda  self,number,base:int(math.ceil(math.log(number,base)))
   __toDecimal= lambda self,binary_num: int(binary_num, 2)

   def __init__(self,image:np.array,kernelSize:int):
    self.model=Model()
    self.__h_params=HGate()
    self.__i_params=IGate()
    self.__x_params=XGate()
    self.__swap_params=SwapGate()
    self.kernelSize=kernelSize


    self.__ireg=QReg(size=8)# grey scale value
    self.__auxreg=QReg(size=1)#helper Qubit

    self.__image=image
    self.__num_rows,self.__num_cols=image.shape
    self.__imageSize=self.__num_rows

    self.r=self.__applyLog(256*kernelSize**2,2)

    self.__rightReg=QReg(size=8)#for addition
    self.__leftargReg=QReg(size=self.r)#for addition
    self.__inreg=QReg(size=self.r)# final value after filter

    if self.__imageSize%2!=0:
      ValueError("Image should be of size 2^n * 2^n ")
      exit

    self.__xreg=QReg(size=int(self.__imageSize/2))#pixel position
    self.__yreg=QReg(size=int(self.__imageSize/2))#pixel position
    self.__encodingpixels=[self.__toBinary(i,self.__imageSize) for i in range(2**self.__imageSize)]

    self.model =self.__initializeQRegister()

   def __initializeQRegister(self):
    model=Model()

    for i in range(8):
      self.__ireg[i]=model.IGate(self.__i_params)[TARGET]
      self.__rightReg[i]=model.IGate(self.__i_params)[TARGET]

    for i in range(int(self.__imageSize/2)):
      self.__xreg[i]=model.HGate(self.__h_params)[TARGET]
      self.__yreg[i]=model.HGate(self.__h_params)[TARGET]
    self.__auxreg=model.IGate(self.__i_params)[TARGET]

    for i in range(self.r):
      self.__inreg[i]=model.IGate(self.__i_params)[TARGET]
      self.__leftargReg[i]=model.IGate(self.__i_params)[TARGET]
    return model


   def authenticate(self):
     classiq.authenticate()

   def encodeToEneqr(self):
    x=int(self.__imageSize/2)
    for m in range(self.__imageSize*self.__imageSize):
      pixel=self.__encodingpixels[m]
      mcx_params=Mcx(num_ctrl_qubits=self.__imageSize,ctrl_state=pixel)

      k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(self.__xreg,self.__yreg),'TARGET_QUBIT':self.__auxreg[0]})
      self.__xreg=k['CTRL_IN'][0:x]
      self.__yreg=k['CTRL_IN'][x:self.__imageSize]
      self.__auxreg[0]=k['TARGET_QUBIT']

      v=m%self.__num_cols
      u=int((m-v)/self.__num_rows)
      pixel_value=self.__toBinary(self.__image[u,v],8)
      for i in range(len(pixel_value)):
          if(pixel_value[i]=='1'):
            cx_params=CXGate(num_ctrl_qubits=1,ctrl_state='1')
            k=self.model.CXGate(cx_params,in_wires={'CTRL':self.__auxreg,'TARGET':self.__ireg[i]})
            self.__auxreg=k['CTRL']
            self.__ireg[i]=k['TARGET']

      k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(self.__xreg,self.__yreg),'TARGET_QUBIT':self.__auxreg[0]})
      self.__xreg=k['CTRL_IN'][0:x]
      self.__yreg=k['CTRL_IN'][x:self.__imageSize]
      self.__auxreg[0]=k['TARGET_QUBIT']



   def encodeToNeqr(self):
    x=int(self.__imageSize/2)
    num_rows,num_cols=self.__image.shape
    for m in range(num_rows*num_cols):
      pixel=self.__encodingpixels[m]
      v=m%num_cols
      u=int((m-v)/num_rows)
      pixel_value=self.__toBinary(self.__image[u,v],8)
      for i in range(len(pixel_value)):
          if(pixel_value[i]=='1'):
            mcx_params=Mcx(num_ctrl_qubits=self.__imageSize,ctrl_state=pixel)
            k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(self.__xreg,self.__yreg),'TARGET_QUBIT':self.__ireg[i]})
            self.__xreg=k['CTRL_IN'][0:x]
            self.__yreg=k['CTRL_IN'][x:self.__imageSize]
            self.__ireg[i]=k['TARGET_QUBIT']

   def cyclicShiftRight(self,sreg:QReg):
    i=1
    x=int(self.__imageSize/2)
    for qubit in range(x,1,-1):
        ctrl_state='1'*(qubit-1)
        num_ctrl_qubits=len(ctrl_state)
        mcx_params=Mcx(num_ctrl_qubits=num_ctrl_qubits,ctrl_state=ctrl_state)
        k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(sreg[i:x]),'TARGET_QUBIT':sreg[i-1]})#changes here
        sreg[i:x]=k['CTRL_IN']
        sreg[i-1]=k['TARGET_QUBIT']
        i=i+1
    k=self.model.XGate(self.__x_params,in_wires={'TARGET':sreg[x-1]})
    sreg[x-1]=k['TARGET']


   def cyclicShiftLeft(self,sreg:QReg):
    x=int(self.__imageSize/2)
    i=x-1
    k=self.model.XGate(self.__x_params,in_wires={'TARGET':sreg[i]})
    sreg[i]=k['TARGET']
    for qubit in range(1,x):
        ctrl_state='1'*qubit
        mcx_params=Mcx(num_ctrl_qubits=qubit,ctrl_state=ctrl_state)
        k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(sreg[i:x]),'TARGET_QUBIT':sreg[x-qubit-1]})#changes here
        sreg[i:x]=k['CTRL_IN']
        sreg[x-qubit-1]=k['TARGET_QUBIT']
        i=i-1


   def swapUp(self): #multiply by 2
    control=int(self.__imageSize/2)
    for i in range(control-1):
        swap_input = QUInt.concat(self.__yreg[control-1],self.__yreg[i])
        k=self.model.SwapGate(self.__swap_params,in_wires={'TARGET':swap_input})
        self.__yreg[control-1]=k['TARGET'][0]
        self.__yreg[i]=k['TARGET'][1]


   def swapdown(self): # x4x3x2x1x0 -> x3x2x1x0x4 #divide by 2
    control=int(self.__imageSize/2)
    for i in range(control):
      if(i==0):
        continue
      swap_input = QUInt.concat(self.__yreg[control-1],self.__yreg[control-1-i])
      k=self.model.SwapGate(self.__swap_params,in_wires={'TARGET':swap_input})
      self.__yreg[control-1]=k['TARGET'][0]
      self.__yreg[control-1-i]=k['TARGET'][1]


   def add(self,x,y):
    self.__controlledCopy(0,0,self.__ireg,self.__rightReg)
    self.__controlledCopy(x,y,self.__inreg,self.__leftargReg,no_of_qubit=self.r)
    add_params=Adder(left_arg={'size':self.r},right_arg={'size':8})
    k=self.model.Adder(add_params,in_wires={'left_arg':self.__leftargReg,'right_arg':self.__rightReg})
    self.__leftargReg=k['left_arg']
    self.__rightReg=k['right_arg']
    self.__controlledCopy(x,y,k['sum'],self.__inreg,no_of_qubit=self.r)

      #upper three step can be reduced if we can reset qubit i.e. qc.reset(self.__inreg) if available


   def __controlledCopy(self,y,x,fromReg:QReg,toReg:QReg,no_of_qubit=8):
        size=int(self.__imageSize/2)
        y=self.__toBinary(y,size)
        x=self.__toBinary(x,size)
        ctrl_state=x+y

        mcx_params=Mcx(num_ctrl_qubits=self.__imageSize,ctrl_state=ctrl_state)
        k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(self.__xreg,self.__yreg),'TARGET_QUBIT':self.__auxreg[0]})
        self.__xreg=k['CTRL_IN'][0:size]
        self.__yreg=k['CTRL_IN'][size:self.__imageSize]
        self.__auxreg[0]=k['TARGET_QUBIT']

        for i in range(no_of_qubit): #perform bitwiseand operation
          ccx_params=CCXGate(num_ctrl_qubits=2,ctrl_state='11')
          k=self.model.CCXGate(ccx_params,in_wires={'CTRL':QUInt.concat(self.__auxreg[0],fromReg[i]),'TARGET':toReg[i]})
          self.__auxreg[0]=k['CTRL'][0]
          fromReg[i]=k['CTRL'][1]
          toReg[i]=k['TARGET']

        k=self.model.Mcx(mcx_params,in_wires={'CTRL_IN':QUInt.concat(self.__xreg,self.__yreg),'TARGET_QUBIT':self.__auxreg[0]})
        self.__xreg=k['CTRL_IN'][0:size]
        self.__yreg=k['CTRL_IN'][size:self.__imageSize]
        self.__auxreg[0]=k['TARGET_QUBIT']

   def apply_filter(self):
    s=int(self.kernelSize/2)
    for y in range(1,self.__imageSize-1):
      for x in range(1,self.__imageSize-1):
        if y<=s-1 or y>=self.__imageSize-s or x<=s-1 or x>=self.__imageSize-s:
          self.__controlledCopy(y,x,self.__ireg,self.__inreg)
        else:
          self.__applyQVXY(x,y)
          # self.__controlledCopy(x,y,self.__vfxReg,self.__inreg)


   def __applyQVXY(self,x,y):
    s=int(self.kernelSize/2)
    for i in range(1,self.kernelSize):
      for j in range(1,self.kernelSize):
        self.add(x,y)
        if j==s:
          for k in range(1,s+1):
            self.cyclicShiftRight(self.__xreg)
        else:
          self.cyclicShiftLeft(self.__xreg)
      if i==s:
        for k in range(1,s+1):
          self.cyclicShiftRight(self.__yreg)
      else:
        self.cyclicShiftLeft(self.__yreg)

   def copy(self):
      self.__controlledCopy(0,0,self.__ireg,self.__inreg)
      self.cyclicShiftRight(self.__yreg)
      self.__controlledCopy(0,0,self.__ireg,self.__inreg)

   def showCircuit(self):
      constraints = Constraints(optimization_parameter=OptimizationParameter.WIDTH)
      self.model.set_outputs({"out": self.__inreg})#since our output is at inreg other not needed
      self.model.sample()
      self.model.constraints=constraints
      self.quantum_program = synthesize(self.model.get_model())
      show(self.quantum_program)

   def arrayToImage(self):
    filterMatrix=self.circuitToImage()
    filteredImage = Image.fromarray(filterMatrix.astype('uint8'))
    return filteredImage.show()

   def circuitToImage(self):# todo call and execute backend and converting binary output to decimal and then executing remaining
    s=int(self.kernelSize/2)
    output=np.array([129,144,156,125,128,129,155,133,119,109,128,153,111,101,151,138])
    output.reshape(self.__imageSize,self.__imageSize)

    for i in range(s,self.__imageSize-s): #adapting classical division to reduce depth of the circuit
      for j in range(s,self.__imageSize-s):
        output[i][j]=output[i][j]/(self.kernelSize**2)

    return output

   def measure(self):
    results = execute(self.quantum_program).result()
