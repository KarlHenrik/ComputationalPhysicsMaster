# for training av validation
import torch
import numpy as np
from vocparseclslabels import PascalVOC
from PIL import Image

# GUI stuff
import PySimpleGUI as sg
from io import BytesIO

# dirs
root_dir = "VOCdevkit/VOC2012/"
img_dir = "VOCdevkit/VOC2012/JPEGImages/"

# getting the correct file names
pVOC = PascalVOC(root_dir)
classes = pVOC.list_image_sets()

# stored values from a stored best model
best_all = torch.load("models/bestMdl.pt")
fnames = np.array(best_all["fnames"])
labels = best_all["labels"].cpu().numpy().astype(int)
preds = best_all["preds"].cpu().numpy()


# function for getting the correct image in png format
def getPNGimage(img_dir, class_name, rank, fnames, preds, classes):
    class_nr = classes.index(class_name)
    sort_inds = preds[:, class_nr].argsort()
    
    fname = fnames[sort_inds][-1 - rank]
    
    filename = img_dir + fname + ".jpg"
    
    im = Image.open(filename)
    bio = BytesIO()
    im.save(bio, format='PNG')
    return bio.getvalue()

# ------------------- GUI -----------------------------
sg.theme("DarkAmber")

img = sg.Image(data = getPNGimage(img_dir, classes[0], 0, fnames, preds, classes), key = "image")

layout = [  [sg.Button("<-----"), sg.Text('0', key = "rank", size=(5, 1), justification="center"), sg.Button("----->"), sg.Button("Close")],
            [sg.Combo(classes, size=(15, 20), default_value = classes[0], enable_events=True, key="combo")],
            [img] ]

# Create the Window
window = sg.Window("Window Title", layout)
rank = 0
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == "Close": # if user closes window or clicks cancel
        break
    if event == "<-----":
        rank -= 1
    if event == "----->":
        rank += 1
    if event == "combo":
        rank = 0
        
    window["rank"].update(str(rank))
    window["image"].update(data = getPNGimage(img_dir, values["combo"], rank, fnames, preds, classes))

window.close()