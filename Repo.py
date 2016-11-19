__author__ = 'Keagan Moo'
    theCountOne = 0
    theCountTwo = 0
    theCountThree = 0
    theCountFour = 0
    theCount = 0
    for i in theMany:
        theSubject = theMany[i]
        thePiece = theSubject.loc['theta']
        thetPiece = Thet[0].loc['theta']
        trueParticle = True
        trueParticleOne = False
        trueParticleTwo = False
        trueParticleThree = False
        trueParticleFour = False
        for j in range(len(thePiece)):
            toMatch = thePiece[j - 1]
            itIs = thetPiece[j - 1]
            if (round(toMatch, 0) != round(itIs, 0)):
                trueParticle = False

            if (round(toMatch, 0) == round(itIs, 0)):
                if (trueParticleOne):
                    if (trueParticleTwo):
                        if(trueParticleThree):
                            trueParticleFour = True

            if (round(toMatch, 0) == round(itIs, 0)):
                if (trueParticleOne):
                    if (trueParticleTwo):
                        trueParticleThree = True

            if (round(toMatch, 0) == round(itIs, 0)):
                if (trueParticleOne):
                    trueParticleTwo = True

            if (round(toMatch, 0) == round(itIs, 0)):
                trueParticleOne = True

        if (trueParticle):
            theCount = theCount + 1
        if (trueParticleOne):
            theCountOne = theCountOne + 1
        if (trueParticleTwo):
            theCountTwo = theCountTwo + 1
        if (trueParticleThree):
            theCountThree = theCountThree + 1
        if (trueParticleFour):
            theCountFour = theCountFour + 1

    print theCount
    print theCountOne
    print theCountTwo
    print theCountThree
    print theCountFour
    theCount = theCount/1000