module metropolis 

export runMetropolis!

using ..harmonicOscillator
using ..slater
using ..gaussian
using ..gaussianSimple
using ..jastrow 
using ..padeJastrow
using ..rbm 
using ..nn

using Statistics
using LinearAlgebra

"""
    runMetropolis!(args...)

Function for running the Metropolis algorithm.

# Arguments
- `system`: The system struct.
- `numMCIterations::Int64`: 
- `stepLength::Float64`: 
- `sampler`:
- `burnIn`:
- `writeToFile`:
- `calculateOnebody`:

# Returns
- `Float`: The local energy of the system. 
"""
function runMetropolis!(
    system, 
    numMcIterations::Int64, 
    stepLength::Float64;
    optimizationIteration::Int64 = 0,
    sampler = "bf", 
    burnIn = 0.001, 
    writeToFile = false, 
    calculateOnebody = false,
    sortInput = false,
    learningRate = 0.0
)
    
    localEnergySum::Float64 = 0.0

    localEnergies::Array{Float64, 1} = zeros(Float64, Int(numMcIterations))

    start = time()

    if sampler == "bf"
        stepFunction = metropolisStepBruteForce!
    elseif sampler == "is"
        stepFunction = metropolisStepImportanceSampling!
    else
        println("Sampler not implemented")
        exit(100)
    end

    if calculateOnebody
        numBins = 6000
        maxLength = 15
        dr = maxLength/numBins
        onebody = zeros(numBins)
    end

    numMoves = 0
    numParticles = system.numParticles

    for i = 1:numMcIterations
        temp = copy(system.particles)
        stepFunction(system, stepLength)

        if sortInput
            sortParticles!(system.particles)
        end

        localEnergy = computeLocalEnergy(system)

        if maximum(abs.(temp-system.particles))!=0
            numMoves += 1
        end

        localEnergies[i] = localEnergy

        if i > ceil(burnIn*numMcIterations)
            for wavefunctionElement in system.wavefunctionElements
                if typeof(wavefunctionElement) != SlaterMatrix
                    wavefunctionElement.variationalParameterGradient = computeParameterGradient(system, wavefunctionElement)
                    wavefunctionElement.localEnergyPsiParameterDerivativeSum += localEnergy*wavefunctionElement.variationalParameterGradient
                    wavefunctionElement.psiParameterDerivativeSum += wavefunctionElement.variationalParameterGradient
                end 
            end
            localEnergySum += localEnergy
                

            if calculateOnebody
                for particle=1:system.numParticles
                    r = sqrt(sum(system.particles[particle,:].^2))
                    if(r <15)
                        onebody[floor(Int, r ÷ dr) + 1] += 1
                    end
                end
            end
        end

        system.iteration += 1

    end


    if calculateOnebody
        filenameOnebody = makeFilenameOnebody(system,stepLength, numMcIterations, sampler, learningRate)
        saveDataToFile(onebody, filenameOnebody)
    end 

    runtime = time() - start

    if writeToFile
        filenameSamples = makeFilenameSamples(system,stepLength, numMcIterations, sampler, learningRate)
        saveDataToFile(localEnergies, filenameSamples)
    end

    samples = numMcIterations - ceil(burnIn*numMcIterations) 

    mcLocalEnergy = localEnergySum/samples

    for wavefunctionElement in system.wavefunctionElements
        if typeof(wavefunctionElement) != SlaterMatrix
            mcLocalEnergyPsiParameterDerivative = wavefunctionElement.localEnergyPsiParameterDerivativeSum/samples
            mcPsiParameterDerivative = wavefunctionElement.psiParameterDerivativeSum/samples
            wavefunctionElement.variationalParameterGradient = 2*(mcLocalEnergyPsiParameterDerivative - mcLocalEnergy*mcPsiParameterDerivative)

            wavefunctionElement.localEnergyPsiParameterDerivativeSum .*=0
            wavefunctionElement.psiParameterDerivativeSum .*=0
        end
    end
    return mcLocalEnergy
end 

"""
    metropolisStepBruteForce!(args...)

Computes the ratio between the current and previous squared wave function value.

# Arguments
- `system`: The system struct.
- `stepLength`:

# Returns
- `Float`: The local energy of the system. 
"""
function metropolisStepBruteForce!(system, stepLength)
    numParticles = system.numParticles 
    numDimensions = system.numDimensions

    # Chooses one coordinate randomly to update.
    coordinateToUpdate::Int64 = rand(1:numDimensions)
    particleToUpdate::Int64 = rand(1:numParticles)

    # Update the coordinate:
    oldPosition = copy(system.particles)
    system.particles[particleToUpdate, coordinateToUpdate] += (rand(Float64) - 0.5)*stepLength

    # Update the slater matrix:
    ratio = 1.0

    for element in system.wavefunctionElements
        updateElement!(system, element, particleToUpdate)
        ratio *= computeRatio(system, 
                            element, 
                            particleToUpdate, 
                            coordinateToUpdate, 
                            oldPosition)
    end

    U = rand(Float64)

    if U < ratio
        if system.slaterInWF
            inverseSlaterMatrixUpdate(system, 
                                    system.wavefunctionElements[1], 
                                    particleToUpdate, 
                                    system.wavefunctionElements[1].R)
        end
    else 
        system.particles[particleToUpdate, coordinateToUpdate] = oldPosition[particleToUpdate, coordinateToUpdate]
        for element in system.wavefunctionElements
            updateElement!(system, element, particleToUpdate)
        end
    end
end

"""
    metropolisStepImportanceSampling!(args...)

Computes the ratio between the current and previous squared wave function value.

# Arguments
- `system`: The system struct.
- `stepLength`:

# Returns
- `Float`: The local energy of the system. 
"""
function metropolisStepImportanceSampling!(system, stepLength)
    numParticles = system.numParticles
    numDimensions = system.numDimensions

    # Chooses one coordinate randomly to update.
    coordinateToUpdate::Int64 = rand(1:numDimensions)
    particleToUpdate::Int64 = rand(1:numParticles)

    # Update the coordinate:
    oldPosition = copy(system.particles)

    D = 0.5

    currentDriftForce = 0.0
    for element in system.wavefunctionElements
        currentDriftForce += computeDriftForce(system, 
                                            element, 
                                            particleToUpdate, 
                                            coordinateToUpdate)
    end

    currentDriftForceFull = zeros(numParticles*numDimensions)
    for element in system.wavefunctionElements
        currentDriftForceFull += 2*computeGradient(system, 
                                            element)
    end

    system.particles[particleToUpdate, coordinateToUpdate] += D*currentDriftForce*stepLength + randn(Float64)*sqrt(stepLength)

    newDriftForce = 0.0
    for element in system.wavefunctionElements
        newDriftForce += computeDriftForce(system, element, particleToUpdate, coordinateToUpdate)
    end

    newDriftForceFull = zeros(numParticles*numDimensions)
    for element in system.wavefunctionElements
        newDriftForceFull += 2*computeGradient(system, 
                                            element)
    end

    greensFunction = computeGreensFunction4(oldPosition,
                                            system.particles,
                                            particleToUpdate,
                                            coordinateToUpdate,
                                            currentDriftForceFull,
                                            newDriftForceFull,
                                            D,
                                            stepLength,
                                            numDimensions)                                   

    ratio = 1.0

    for element in system.wavefunctionElements
        updateElement!(system, element, particleToUpdate)
        ratio *= computeRatio(system, 
                            element,
                            particleToUpdate,
                            coordinateToUpdate,
                            oldPosition)
    end

    U = rand(Float64)



    if U < greensFunction*ratio
        if system.slaterInWF
            inverseSlaterMatrixUpdate(system, system.wavefunctionElements[1], particleToUpdate, system.wavefunctionElements[1].R)
        end
    else 
        system.particles[particleToUpdate, coordinateToUpdate] = oldPosition[particleToUpdate, coordinateToUpdate]
        for element in system.wavefunctionElements
            updateElement!(system, element, particleToUpdate)
        end
    end
end

"""
    computeGreensFunction(args...)

Computes the ratio between the current and previous squared wave function value.

# Arguments
- `system`: The system struct.
- `stepLength`:

# Returns
- `Float`: The local energy of the system. 
"""
# function computeGreensFunction(oldPosition, 
#                             newPosition, 
#                             particleToUpdate,        
#                             coordinateToUpdate, 
#                             oldDriftForce, 
#                             newDriftForce, 
#                             D,
#                             stepLength)

#     greensFunctionArgument = (oldPosition[particleToUpdate, coordinateToUpdate] +
#                                 - newPosition[particleToUpdate, coordinateToUpdate] +
#                                 - D*stepLength*newDriftForce)^2 +
#                                 - (newPosition[particleToUpdate, coordinateToUpdate] + 
#                                 - oldPosition[particleToUpdate, coordinateToUpdate] +
#                                 - D*stepLength*oldDriftForce)^2

#     greensFunctionArgument /= (4.0*D*stepLength)
#     greensFunction = exp(-greensFunctionArgument)
#     return greensFunction
# end

function computeGreensFunction2(oldPosition, 
                            newPosition, 
                            oldDriftForce, 
                            newDriftForce, 
                            particleToUpdate,
                            D,
                            stepLength,
                            numDims)

    greensFunction1 = exp(-dot(oldPosition[particleToUpdate,:] - newPosition[particleToUpdate,:] - D*stepLength*newDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims], oldPosition[particleToUpdate,:] - newPosition[particleToUpdate,:] - D*stepLength*newDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims])/(4*D*stepLength))
    greensFunction2 = exp(-dot(newPosition[particleToUpdate,:] - oldPosition[particleToUpdate,:] - D*stepLength*oldDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims], newPosition[particleToUpdate,:] - oldPosition[particleToUpdate,:] - D*stepLength*oldDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims])/(4*D*stepLength))
    
    return greensFunction1/greensFunction2

end

function computeGreensFunction3(oldPosition, 
                            newPosition, 
                            oldDriftForce, 
                            newDriftForce, 
                            particleToUpdate,
                            D,
                            stepLength)

    greensFunction1 = exp(-dot(reshape(oldPosition',1,:)' - reshape(newPosition',1,:)' - D*stepLength*newDriftForce, reshape(oldPosition',1,:)' - reshape(newPosition',1,:)' - D*stepLength*newDriftForce)/(4*D*stepLength))
    greensFunction2 = exp(-dot(reshape(newPosition',1,:)' - reshape(oldPosition',1,:)' - D*stepLength*oldDriftForce, reshape(newPosition',1,:)' - reshape(oldPosition',1,:)' - D*stepLength*oldDriftForce)/(4*D*stepLength))

    return greensFunction1/greensFunction2

end

function computeGreensFunction4(oldPosition, 
                            newPosition, 
                            particleToUpdate,
                            coordinateToUpdate,
                            oldDriftForce, 
                            newDriftForce, 
                            D,
                            stepLength,
                            numDims)

    arg = 0
    oldDriftForce = oldDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims]
    newDriftForce = newDriftForce[(particleToUpdate-1)*numDims + 1: (particleToUpdate-1)*numDims+numDims]
    for i=1:numDims
        arg += 0.5*(oldDriftForce[i] + newDriftForce[i])*(D*stepLength*0.5*(oldDriftForce[i] - newDriftForce[i]) - newPosition[particleToUpdate, i] + oldPosition[particleToUpdate, i])
    end
    return exp(arg)

end

function computeGreensFunction5(oldPosition, 
    newPosition, 
    oldDriftForce, 
    newDriftForce, 
    particleToUpdate,
    coordinateToUpdate,
    D,
    stepLength)
    return exp(0.5*(oldDriftForce - newDriftForce)*(newPosition[particleToUpdate, coordinateToUpdate] -oldPosition[particleToUpdate, coordinateToUpdate]) ) + 1
end

function saveDataToFile(data, filename::String)
    return
end

function makeFilenameSamples(system, steplength, numMCsteps, sampler, learningRate)
    wavefunctionCombination = "wf_"
    wavefunctionElementsInfo = "_elementinfo_"
    for element in system.wavefunctionElements
        elementinfo = wavefunctionName(element) 
        wavefunctionCombination  = wavefunctionCombination * elementinfo[2] * "_"
        wavefunctionElementsInfo = wavefunctionElementsInfo * elementinfo[1] * "_"
    end

    if system.interacting == true
        folder = "Interacting"
    elseif system.interacting == false
        folder = "Non_Interacting"
    end
    filename = "Data/"* system.hamiltonian * "/MC/" * folder * "/" * wavefunctionCombination * "sysInfo_" * sampler  *"_omega_" * string(system.omega) * "_sl_" * string(steplength)* "_mcSteps_"* string(numMCsteps) * "_lr_" * string(learningRate) * "_numD_" * string(system.numDimensions) * "_numP_" * string(system.numParticles) * wavefunctionElementsInfo *".txt"
    return filename
end

function makeFilenameOnebody(system, steplength, numMCsteps, sampler, learningRate)
    wavefunctionCombination = "wf_"
    wavefunctionElementsInfo = "_elementinfo_"
    for element in system.wavefunctionElements
        elementinfo = wavefunctionName(element) 
        wavefunctionCombination  = wavefunctionCombination * elementinfo[2] * "_"
        wavefunctionElementsInfo = wavefunctionElementsInfo * elementinfo[1] * "_"
    end

    if system.interacting == true
        folder = "Interacting"
    elseif system.interacting == false
        folder = "Non_Interacting"
    end
    filename = "Data/"* system.hamiltonian * "/Onebody/" * folder * "/" * wavefunctionCombination * "sysInfo_" * sampler  *"_omega_" * string(system.omega) * "_sl_" * string(steplength)* "_mcSteps_"* string(numMCsteps) * "_lr_" * string(learningRate) * "_numD_" * string(system.numDimensions) * "_numP_" * string(system.numParticles) * wavefunctionElementsInfo *".txt"
    return filename
end

function wavefunctionName(element::SlaterMatrix)
    return ["slater_none", "slater"]
end

function wavefunctionName(element::Gaussian)
    return ["gaussian_none", "gaussian"]
end

function wavefunctionName(element::Jastrow)
    return ["jastrow_none", "jastrow"]
end

function wavefunctionName(element::PadeJastrow)
    return ["padeJastrow_none", "padeJastrow"]
end

function wavefunctionName(element::GaussianSimple)
    return ["gaussianSimple_none", "gaussianSimple"]
end

function wavefunctionName(element::RBM)
    return [("rbm_numhidden_" * string(size(element.h)[1])), "rbm"]
end

function wavefunctionName(element::NN)
    return [("nn_nh1_" * string(size(element.a[1])[1]) * "_nh2_" * string(size(element.a[2])[1])) * "_af_" * string(element.activationFunction), "nn"]
end


function sortParticles!(particles)
    distance = vec(sqrt.(sum(particles.^2,dims=2)))
    perm = sortperm(distance)
    particles[:,:] = particles[perm, :]
end

end