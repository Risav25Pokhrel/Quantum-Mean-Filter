<h1>QUANTUM MEAN FILTER USING CLASSIQ</h1>

<h2>Running on Device</h2>
!pip install -U classiq

<h3>Calling Method</h3>

sample_image=ImageToArray().getImage("...path",n=2)
meanfilter=ImageToCircuit(sample_image,2)
meanfilter.authenticate() #should authenticate if your device is not authenticate otherwise ignore this step
meanfilter.encodeToEneqr()
meanfilter.apply_filter()
meanfilter.showCircuit()


<h2>NEQR and ENEQR of 2 x 2 Image </h2>
![image](https://github.com/Risav25Pokhrel/Quantum-Mean-Filter/assets/103576193/f48fa04d-4783-46e9-98cd-42320b662969)

<h2>Cyclic shift and Copy Results</h2>
![image](https://github.com/Risav25Pokhrel/Quantum-Mean-Filter/assets/103576193/e3886f27-3ed7-45d2-bfa6-3882196e26e8)


<h2>Final Image of Mean Filter</h2>
![image](https://github.com/Risav25Pokhrel/Quantum-Mean-Filter/assets/103576193/74c5087e-983c-4d0b-abdf-5074e0e0a987)



TODO:
 Applying mean filter for GQIR image representation for any size image
 Implementation of the Quantum Non-Modular Floating Point Division
 Normalization and Filtering of output 
 
