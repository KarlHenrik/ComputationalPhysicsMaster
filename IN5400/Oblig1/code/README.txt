There are three runnable python files. They are run with no command line parameters.

trainModel.py -> trains a Multi-label prediction network, saves the best one to the models folder,
		 and shows some plots.

demonstrator.py -> loads the best model from trainModel.py and lets you see the images ordered
		   according to their prediction, for all classes, in a GUI.
	           You will need to install pysimplegui.

task2.py -> creates two models. where one has some custom convolutions and batchnorms incorperationg
	    some normalization. predicts on some data and prints how close their predictions were.
	    We gave up on training the custom model.

The other python files are helper functions given by the instructors, except for vocTrainEval.py,
which includes the classes and functions we wrote to be able to train the model in trainModel.py.

The VOC data must be located in the same folder as the .py files. We assume there is a folder "VOCdevkit"
which includes a folder VOC2012, which again includes folders for the JPEG images and labels.