#ifndef HcalTrigPrimRecHitProducer_h
#define HcalTrigPrimRecHitProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalCoderFactory.h"
#include "SimCalorimetry/HcalSimAlgos/interface/HcalTriggerPrimitiveAlgo.h"

using namespace cms;

class HcalTrigPrimRecHitProducer : public edm::EDProducer
{
public:

  explicit HcalTrigPrimRecHitProducer(const edm::ParameterSet& ps);
  virtual ~HcalTrigPrimRecHitProducer() {}

  /**Produces the EDM products,*/
  virtual void produce(edm::Event& e, const edm::EventSetup& c);

private:

  HcalCoderFactory theCoderFactory;
  HcalTriggerPrimitiveAlgo theAlgo;
};

#endif

