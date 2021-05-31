#include "analyzePythia.h"

void analyzePythia(TString inputFile, TString outputFile) {

    //Loading everything
    AliRunLoader* inputRun = AliRunLoader::Open(inputFile);
    int numEvents = inputRun->GetNumberOfEvents();
    inputRun->LoadKinematics();
    inputRun->LoadHeader();

    //Loop through all of the events in the run
    for(int event = 0;  event < numEvents; event++) {

        //Load the stack (contains all of the particles in the event)
        inputRun->GetEvent(event);
        AliStack *theStack = inputRun->Stack();

        int numPrimary = theStack->GetNprimary();
        int numTracks = theStack->GetNtrack();

        std::cout << numPrimary << " is num primary and " << numTracks << " is num tracks\n";

        for(int part = 0; part < numPrimary; part++) {

            //Load the particle and do something with it
            TParticle *particle = theStack->Particle(part);

            // If particle is a charged hadron
            if(std::find(std::begin(charged_pdg_list), std::end(charged_pdg_list), particle->GetPdgCode()) != std::end(charged_pdg_list)) {
                float distPoint[4];
                distPoint[0] = particle->Pt();
                distPoint[1] = particle->Phi();
                distPoint[2] = particle->Eta();
                distPoint[3] = particle->GetPdgCode();
            }
        }
    }
}