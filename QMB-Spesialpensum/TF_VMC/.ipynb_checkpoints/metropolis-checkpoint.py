from WaveFunctions import *
import tensorflow as tf
import numpy as np

def metro_util(config):
    # Creating the correct step function
    metro_name = config["metropolis"]["name"]
    if metro_name == "metro":
        step = create_metro_step(config)
    elif metro_name == "importance":
        step = create_importance_step(config)
    else:
        raise NotImplementedError
        
    return step

"""
def create_metro_step(config):
    num = config["particles"]
    dims = config["dimensions"]
    step_length = config["metropolis"]["step_length"]
    
    def metro_step(particles, wf):
        p1 = np.random.randint(num)
        oldPos = particles.getPos(p1)

        adj_dir = np.random.randint(dims)
        adj_sign = (np.random.randint(2) - 0.5) * 2
        newPos = oldPos + 0
        newPos[adj_dir] += step_length * adj_sign

        particles.setPos(p1, newPos)
        ratio = wf.ratio(particles, p1, oldPos)
        if (np.random.random() < np.power(ratio, 2)):
            return True
        else:
            particles.setPos(p1, oldPos)
            return False
    return metro_step
"""

def create_importance_step(config):
    num = tf.constant(config["particles"])
    dims = tf.constant(config["dimensions"])
    step_length = tf.constant(config["metropolis"]["step_length"], dtype = tf.float64)
    
    one_once = np.zeros((num, 1))
    one_once[0] = 1
    one_once = tf.constant(one_once)

    @tf.function
    def importance_step(old_positions, wf):
        chosen_particle_mask = tf.random.shuffle(one_once)
        p1 = tf.math.argmax(chosen_particle_mask)[0]
        oldQF = wf.QF(old_positions, p1)
        
        move = 0.5 * oldQF * step_length + tf.random.normal(shape = (dims,), dtype = tf.float64) * tf.math.sqrt(step_length)
        new_positions = old_positions + chosen_particle_mask * move
        
        ratio = wf.ratio(new_positions, p1, old_positions)
        newQF = wf.QF(new_positions, p1)

        greensFuncRatio = tf.math.reduce_sum( (oldQF + newQF) * (old_positions[p1] - new_positions[p1] + step_length * (oldQF - newQF)) )
        greensFuncRatio = tf.math.exp(0.5 * greensFuncRatio)
        if (tf.random.uniform(shape = (), dtype = tf.float64) < greensFuncRatio * tf.math.pow(ratio, 2)):
            return new_positions
        else:
            return old_positions
    return importance_step