// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdIsrcdIWCLRootDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "WCSimRecoObjectTable.hh"
#include "WCSimRecoDigit.hh"
#include "SandBoxPMTcoverage.hh"
#include "WCSimTrueLight.hh"
#include "WCSimTruePart.hh"
#include "WCSimTrueCapture.hh"
#include "WCLTreeReader.hh"
#include "WCLTreeWriter.hh"
#include "WChRecoLite.hh"
#include "HighEReco.hh"
#include "LowEReco.hh"
#include "LikelihoodGenerator.hh"

// Header files passed via #pragma extra_include

namespace ROOT {

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoObjectTable*)
   {
      ::WCSimRecoObjectTable *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoObjectTable >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoObjectTable", ::WCSimRecoObjectTable::Class_Version(), "WCSimRecoObjectTable.hh", 6,
                  typeid(::WCSimRecoObjectTable), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRecoObjectTable::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoObjectTable) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoObjectTable*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoObjectTable*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCSimRecoDigit(void *p);
   static void deleteArray_WCSimRecoDigit(void *p);
   static void destruct_WCSimRecoDigit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRecoDigit*)
   {
      ::WCSimRecoDigit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRecoDigit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRecoDigit", ::WCSimRecoDigit::Class_Version(), "WCSimRecoDigit.hh", 6,
                  typeid(::WCSimRecoDigit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRecoDigit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRecoDigit) );
      instance.SetDelete(&delete_WCSimRecoDigit);
      instance.SetDeleteArray(&deleteArray_WCSimRecoDigit);
      instance.SetDestructor(&destruct_WCSimRecoDigit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRecoDigit*)
   {
      return GenerateInitInstanceLocal((::WCSimRecoDigit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WChRecoLite*)
   {
      ::WChRecoLite *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WChRecoLite >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WChRecoLite", ::WChRecoLite::Class_Version(), "WChRecoLite.hh", 34,
                  typeid(::WChRecoLite), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WChRecoLite::Dictionary, isa_proxy, 4,
                  sizeof(::WChRecoLite) );
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WChRecoLite*)
   {
      return GenerateInitInstanceLocal((::WChRecoLite*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WChRecoLite*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCLTreeReader(void *p);
   static void deleteArray_WCLTreeReader(void *p);
   static void destruct_WCLTreeReader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCLTreeReader*)
   {
      ::WCLTreeReader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCLTreeReader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCLTreeReader", ::WCLTreeReader::Class_Version(), "WCLTreeReader.hh", 15,
                  typeid(::WCLTreeReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCLTreeReader::Dictionary, isa_proxy, 4,
                  sizeof(::WCLTreeReader) );
      instance.SetDelete(&delete_WCLTreeReader);
      instance.SetDeleteArray(&deleteArray_WCLTreeReader);
      instance.SetDestructor(&destruct_WCLTreeReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCLTreeReader*)
   {
      return GenerateInitInstanceLocal((::WCLTreeReader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCLTreeReader*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCLTreeWriter(void *p);
   static void deleteArray_WCLTreeWriter(void *p);
   static void destruct_WCLTreeWriter(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCLTreeWriter*)
   {
      ::WCLTreeWriter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCLTreeWriter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCLTreeWriter", ::WCLTreeWriter::Class_Version(), "WCLTreeWriter.hh", 17,
                  typeid(::WCLTreeWriter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCLTreeWriter::Dictionary, isa_proxy, 4,
                  sizeof(::WCLTreeWriter) );
      instance.SetDelete(&delete_WCLTreeWriter);
      instance.SetDeleteArray(&deleteArray_WCLTreeWriter);
      instance.SetDestructor(&destruct_WCLTreeWriter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCLTreeWriter*)
   {
      return GenerateInitInstanceLocal((::WCLTreeWriter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCLTreeWriter*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SandBoxPMTcoverage(void *p = 0);
   static void *newArray_SandBoxPMTcoverage(Long_t size, void *p);
   static void delete_SandBoxPMTcoverage(void *p);
   static void deleteArray_SandBoxPMTcoverage(void *p);
   static void destruct_SandBoxPMTcoverage(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SandBoxPMTcoverage*)
   {
      ::SandBoxPMTcoverage *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SandBoxPMTcoverage >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SandBoxPMTcoverage", ::SandBoxPMTcoverage::Class_Version(), "SandBoxPMTcoverage.hh", 6,
                  typeid(::SandBoxPMTcoverage), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SandBoxPMTcoverage::Dictionary, isa_proxy, 4,
                  sizeof(::SandBoxPMTcoverage) );
      instance.SetNew(&new_SandBoxPMTcoverage);
      instance.SetNewArray(&newArray_SandBoxPMTcoverage);
      instance.SetDelete(&delete_SandBoxPMTcoverage);
      instance.SetDeleteArray(&deleteArray_SandBoxPMTcoverage);
      instance.SetDestructor(&destruct_SandBoxPMTcoverage);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SandBoxPMTcoverage*)
   {
      return GenerateInitInstanceLocal((::SandBoxPMTcoverage*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SandBoxPMTcoverage*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCSimTrueLight(void *p);
   static void deleteArray_WCSimTrueLight(void *p);
   static void destruct_WCSimTrueLight(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTrueLight*)
   {
      ::WCSimTrueLight *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTrueLight >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTrueLight", ::WCSimTrueLight::Class_Version(), "WCSimTrueLight.hh", 7,
                  typeid(::WCSimTrueLight), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimTrueLight::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTrueLight) );
      instance.SetDelete(&delete_WCSimTrueLight);
      instance.SetDeleteArray(&deleteArray_WCSimTrueLight);
      instance.SetDestructor(&destruct_WCSimTrueLight);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTrueLight*)
   {
      return GenerateInitInstanceLocal((::WCSimTrueLight*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTrueLight*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCSimTruePart(void *p);
   static void deleteArray_WCSimTruePart(void *p);
   static void destruct_WCSimTruePart(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTruePart*)
   {
      ::WCSimTruePart *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTruePart >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTruePart", ::WCSimTruePart::Class_Version(), "WCSimTruePart.hh", 7,
                  typeid(::WCSimTruePart), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimTruePart::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTruePart) );
      instance.SetDelete(&delete_WCSimTruePart);
      instance.SetDeleteArray(&deleteArray_WCSimTruePart);
      instance.SetDestructor(&destruct_WCSimTruePart);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTruePart*)
   {
      return GenerateInitInstanceLocal((::WCSimTruePart*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTruePart*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WCSimTrueCapture(void *p);
   static void deleteArray_WCSimTrueCapture(void *p);
   static void destruct_WCSimTrueCapture(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimTrueCapture*)
   {
      ::WCSimTrueCapture *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimTrueCapture >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimTrueCapture", ::WCSimTrueCapture::Class_Version(), "WCSimTrueCapture.hh", 7,
                  typeid(::WCSimTrueCapture), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimTrueCapture::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimTrueCapture) );
      instance.SetDelete(&delete_WCSimTrueCapture);
      instance.SetDeleteArray(&deleteArray_WCSimTrueCapture);
      instance.SetDestructor(&destruct_WCSimTrueCapture);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimTrueCapture*)
   {
      return GenerateInitInstanceLocal((::WCSimTrueCapture*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimTrueCapture*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_HighEReco(void *p = 0);
   static void *newArray_HighEReco(Long_t size, void *p);
   static void delete_HighEReco(void *p);
   static void deleteArray_HighEReco(void *p);
   static void destruct_HighEReco(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HighEReco*)
   {
      ::HighEReco *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HighEReco >(0);
      static ::ROOT::TGenericClassInfo 
         instance("HighEReco", ::HighEReco::Class_Version(), "HighEReco.hh", 7,
                  typeid(::HighEReco), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HighEReco::Dictionary, isa_proxy, 4,
                  sizeof(::HighEReco) );
      instance.SetNew(&new_HighEReco);
      instance.SetNewArray(&newArray_HighEReco);
      instance.SetDelete(&delete_HighEReco);
      instance.SetDeleteArray(&deleteArray_HighEReco);
      instance.SetDestructor(&destruct_HighEReco);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HighEReco*)
   {
      return GenerateInitInstanceLocal((::HighEReco*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::HighEReco*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LowEReco(void *p = 0);
   static void *newArray_LowEReco(Long_t size, void *p);
   static void delete_LowEReco(void *p);
   static void deleteArray_LowEReco(void *p);
   static void destruct_LowEReco(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LowEReco*)
   {
      ::LowEReco *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LowEReco >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LowEReco", ::LowEReco::Class_Version(), "LowEReco.hh", 6,
                  typeid(::LowEReco), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LowEReco::Dictionary, isa_proxy, 4,
                  sizeof(::LowEReco) );
      instance.SetNew(&new_LowEReco);
      instance.SetNewArray(&newArray_LowEReco);
      instance.SetDelete(&delete_LowEReco);
      instance.SetDeleteArray(&deleteArray_LowEReco);
      instance.SetDestructor(&destruct_LowEReco);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LowEReco*)
   {
      return GenerateInitInstanceLocal((::LowEReco*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LowEReco*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_LikelihoodGenerator(void *p = 0);
   static void *newArray_LikelihoodGenerator(Long_t size, void *p);
   static void delete_LikelihoodGenerator(void *p);
   static void deleteArray_LikelihoodGenerator(void *p);
   static void destruct_LikelihoodGenerator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::LikelihoodGenerator*)
   {
      ::LikelihoodGenerator *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::LikelihoodGenerator >(0);
      static ::ROOT::TGenericClassInfo 
         instance("LikelihoodGenerator", ::LikelihoodGenerator::Class_Version(), "LikelihoodGenerator.hh", 23,
                  typeid(::LikelihoodGenerator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::LikelihoodGenerator::Dictionary, isa_proxy, 4,
                  sizeof(::LikelihoodGenerator) );
      instance.SetNew(&new_LikelihoodGenerator);
      instance.SetNewArray(&newArray_LikelihoodGenerator);
      instance.SetDelete(&delete_LikelihoodGenerator);
      instance.SetDeleteArray(&deleteArray_LikelihoodGenerator);
      instance.SetDestructor(&destruct_LikelihoodGenerator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::LikelihoodGenerator*)
   {
      return GenerateInitInstanceLocal((::LikelihoodGenerator*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::LikelihoodGenerator*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr WCSimRecoObjectTable::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoObjectTable::Class_Name()
{
   return "WCSimRecoObjectTable";
}

//______________________________________________________________________________
const char *WCSimRecoObjectTable::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoObjectTable::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRecoObjectTable::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoObjectTable::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoObjectTable*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRecoDigit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRecoDigit::Class_Name()
{
   return "WCSimRecoDigit";
}

//______________________________________________________________________________
const char *WCSimRecoDigit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRecoDigit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRecoDigit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRecoDigit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRecoDigit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WChRecoLite::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WChRecoLite::Class_Name()
{
   return "WChRecoLite";
}

//______________________________________________________________________________
const char *WChRecoLite::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WChRecoLite*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WChRecoLite::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WChRecoLite*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WChRecoLite::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WChRecoLite*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WChRecoLite::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WChRecoLite*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCLTreeReader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCLTreeReader::Class_Name()
{
   return "WCLTreeReader";
}

//______________________________________________________________________________
const char *WCLTreeReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeReader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCLTreeReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeReader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCLTreeReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeReader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCLTreeReader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeReader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCLTreeWriter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCLTreeWriter::Class_Name()
{
   return "WCLTreeWriter";
}

//______________________________________________________________________________
const char *WCLTreeWriter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeWriter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCLTreeWriter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeWriter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCLTreeWriter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeWriter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCLTreeWriter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCLTreeWriter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SandBoxPMTcoverage::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SandBoxPMTcoverage::Class_Name()
{
   return "SandBoxPMTcoverage";
}

//______________________________________________________________________________
const char *SandBoxPMTcoverage::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SandBoxPMTcoverage*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SandBoxPMTcoverage::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SandBoxPMTcoverage*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SandBoxPMTcoverage::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SandBoxPMTcoverage*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SandBoxPMTcoverage::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SandBoxPMTcoverage*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimTrueLight::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTrueLight::Class_Name()
{
   return "WCSimTrueLight";
}

//______________________________________________________________________________
const char *WCSimTrueLight::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueLight*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTrueLight::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueLight*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimTrueLight::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueLight*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTrueLight::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueLight*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimTruePart::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTruePart::Class_Name()
{
   return "WCSimTruePart";
}

//______________________________________________________________________________
const char *WCSimTruePart::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTruePart*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTruePart::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTruePart*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimTruePart::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTruePart*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTruePart::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTruePart*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimTrueCapture::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimTrueCapture::Class_Name()
{
   return "WCSimTrueCapture";
}

//______________________________________________________________________________
const char *WCSimTrueCapture::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueCapture*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimTrueCapture::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueCapture*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimTrueCapture::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueCapture*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimTrueCapture::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimTrueCapture*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HighEReco::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *HighEReco::Class_Name()
{
   return "HighEReco";
}

//______________________________________________________________________________
const char *HighEReco::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HighEReco*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int HighEReco::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HighEReco*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HighEReco::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HighEReco*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HighEReco::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HighEReco*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LowEReco::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LowEReco::Class_Name()
{
   return "LowEReco";
}

//______________________________________________________________________________
const char *LowEReco::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LowEReco*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LowEReco::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LowEReco*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LowEReco::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LowEReco*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LowEReco::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LowEReco*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr LikelihoodGenerator::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *LikelihoodGenerator::Class_Name()
{
   return "LikelihoodGenerator";
}

//______________________________________________________________________________
const char *LikelihoodGenerator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodGenerator*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int LikelihoodGenerator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodGenerator*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *LikelihoodGenerator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodGenerator*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *LikelihoodGenerator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::LikelihoodGenerator*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void WCSimRecoObjectTable::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoObjectTable.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoObjectTable::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoObjectTable::Class(),this);
   }
}

namespace ROOT {
} // end of namespace ROOT for class ::WCSimRecoObjectTable

//______________________________________________________________________________
void WCSimRecoDigit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRecoDigit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRecoDigit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRecoDigit::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimRecoDigit(void *p) {
      delete ((::WCSimRecoDigit*)p);
   }
   static void deleteArray_WCSimRecoDigit(void *p) {
      delete [] ((::WCSimRecoDigit*)p);
   }
   static void destruct_WCSimRecoDigit(void *p) {
      typedef ::WCSimRecoDigit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRecoDigit

//______________________________________________________________________________
void WChRecoLite::Streamer(TBuffer &R__b)
{
   // Stream an object of class WChRecoLite.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WChRecoLite::Class(),this);
   } else {
      R__b.WriteClassBuffer(WChRecoLite::Class(),this);
   }
}

namespace ROOT {
} // end of namespace ROOT for class ::WChRecoLite

//______________________________________________________________________________
void WCLTreeReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCLTreeReader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCLTreeReader::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCLTreeReader::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCLTreeReader(void *p) {
      delete ((::WCLTreeReader*)p);
   }
   static void deleteArray_WCLTreeReader(void *p) {
      delete [] ((::WCLTreeReader*)p);
   }
   static void destruct_WCLTreeReader(void *p) {
      typedef ::WCLTreeReader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCLTreeReader

//______________________________________________________________________________
void WCLTreeWriter::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCLTreeWriter.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCLTreeWriter::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCLTreeWriter::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCLTreeWriter(void *p) {
      delete ((::WCLTreeWriter*)p);
   }
   static void deleteArray_WCLTreeWriter(void *p) {
      delete [] ((::WCLTreeWriter*)p);
   }
   static void destruct_WCLTreeWriter(void *p) {
      typedef ::WCLTreeWriter current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCLTreeWriter

//______________________________________________________________________________
void SandBoxPMTcoverage::Streamer(TBuffer &R__b)
{
   // Stream an object of class SandBoxPMTcoverage.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SandBoxPMTcoverage::Class(),this);
   } else {
      R__b.WriteClassBuffer(SandBoxPMTcoverage::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SandBoxPMTcoverage(void *p) {
      return  p ? new(p) ::SandBoxPMTcoverage : new ::SandBoxPMTcoverage;
   }
   static void *newArray_SandBoxPMTcoverage(Long_t nElements, void *p) {
      return p ? new(p) ::SandBoxPMTcoverage[nElements] : new ::SandBoxPMTcoverage[nElements];
   }
   // Wrapper around operator delete
   static void delete_SandBoxPMTcoverage(void *p) {
      delete ((::SandBoxPMTcoverage*)p);
   }
   static void deleteArray_SandBoxPMTcoverage(void *p) {
      delete [] ((::SandBoxPMTcoverage*)p);
   }
   static void destruct_SandBoxPMTcoverage(void *p) {
      typedef ::SandBoxPMTcoverage current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SandBoxPMTcoverage

//______________________________________________________________________________
void WCSimTrueLight::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTrueLight.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTrueLight::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTrueLight::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTrueLight(void *p) {
      delete ((::WCSimTrueLight*)p);
   }
   static void deleteArray_WCSimTrueLight(void *p) {
      delete [] ((::WCSimTrueLight*)p);
   }
   static void destruct_WCSimTrueLight(void *p) {
      typedef ::WCSimTrueLight current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTrueLight

//______________________________________________________________________________
void WCSimTruePart::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTruePart.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTruePart::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTruePart::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTruePart(void *p) {
      delete ((::WCSimTruePart*)p);
   }
   static void deleteArray_WCSimTruePart(void *p) {
      delete [] ((::WCSimTruePart*)p);
   }
   static void destruct_WCSimTruePart(void *p) {
      typedef ::WCSimTruePart current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTruePart

//______________________________________________________________________________
void WCSimTrueCapture::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimTrueCapture.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimTrueCapture::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimTrueCapture::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WCSimTrueCapture(void *p) {
      delete ((::WCSimTrueCapture*)p);
   }
   static void deleteArray_WCSimTrueCapture(void *p) {
      delete [] ((::WCSimTrueCapture*)p);
   }
   static void destruct_WCSimTrueCapture(void *p) {
      typedef ::WCSimTrueCapture current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimTrueCapture

//______________________________________________________________________________
void HighEReco::Streamer(TBuffer &R__b)
{
   // Stream an object of class HighEReco.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(HighEReco::Class(),this);
   } else {
      R__b.WriteClassBuffer(HighEReco::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HighEReco(void *p) {
      return  p ? new(p) ::HighEReco : new ::HighEReco;
   }
   static void *newArray_HighEReco(Long_t nElements, void *p) {
      return p ? new(p) ::HighEReco[nElements] : new ::HighEReco[nElements];
   }
   // Wrapper around operator delete
   static void delete_HighEReco(void *p) {
      delete ((::HighEReco*)p);
   }
   static void deleteArray_HighEReco(void *p) {
      delete [] ((::HighEReco*)p);
   }
   static void destruct_HighEReco(void *p) {
      typedef ::HighEReco current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HighEReco

//______________________________________________________________________________
void LowEReco::Streamer(TBuffer &R__b)
{
   // Stream an object of class LowEReco.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LowEReco::Class(),this);
   } else {
      R__b.WriteClassBuffer(LowEReco::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LowEReco(void *p) {
      return  p ? new(p) ::LowEReco : new ::LowEReco;
   }
   static void *newArray_LowEReco(Long_t nElements, void *p) {
      return p ? new(p) ::LowEReco[nElements] : new ::LowEReco[nElements];
   }
   // Wrapper around operator delete
   static void delete_LowEReco(void *p) {
      delete ((::LowEReco*)p);
   }
   static void deleteArray_LowEReco(void *p) {
      delete [] ((::LowEReco*)p);
   }
   static void destruct_LowEReco(void *p) {
      typedef ::LowEReco current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LowEReco

//______________________________________________________________________________
void LikelihoodGenerator::Streamer(TBuffer &R__b)
{
   // Stream an object of class LikelihoodGenerator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(LikelihoodGenerator::Class(),this);
   } else {
      R__b.WriteClassBuffer(LikelihoodGenerator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_LikelihoodGenerator(void *p) {
      return  p ? new(p) ::LikelihoodGenerator : new ::LikelihoodGenerator;
   }
   static void *newArray_LikelihoodGenerator(Long_t nElements, void *p) {
      return p ? new(p) ::LikelihoodGenerator[nElements] : new ::LikelihoodGenerator[nElements];
   }
   // Wrapper around operator delete
   static void delete_LikelihoodGenerator(void *p) {
      delete ((::LikelihoodGenerator*)p);
   }
   static void deleteArray_LikelihoodGenerator(void *p) {
      delete [] ((::LikelihoodGenerator*)p);
   }
   static void destruct_LikelihoodGenerator(void *p) {
      typedef ::LikelihoodGenerator current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::LikelihoodGenerator

namespace {
  void TriggerDictionaryInitialization_WCLRootDict_Impl() {
    static const char* headers[] = {
"WCSimRecoObjectTable.hh",
"WCSimRecoDigit.hh",
"SandBoxPMTcoverage.hh",
"WCSimTrueLight.hh",
"WCSimTruePart.hh",
"WCSimTrueCapture.hh",
"WCLTreeReader.hh",
"WCLTreeWriter.hh",
"WChRecoLite.hh",
"HighEReco.hh",
"LowEReco.hh",
"LikelihoodGenerator.hh",
0
    };
    static const char* includePaths[] = {
"./include",
"/include",
"/exports/applications/apps/SL7/root/6.06.02/include",
"/exports/eddie3_apps_local/apps/SL7/root/6.06.02/include",
"/exports/csce/datastore/ph/groups/PPE/titus/ts-WChRecoSandBox/scripts/nicolas/optimised/source4/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "WCLRootDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$WCSimRecoObjectTable.hh")))  WCSimRecoObjectTable;
class __attribute__((annotate("$clingAutoload$WCSimRecoDigit.hh")))  WCSimRecoDigit;
class __attribute__((annotate("$clingAutoload$WChRecoLite.hh")))  WChRecoLite;
class __attribute__((annotate("$clingAutoload$WCLTreeReader.hh")))  WCLTreeReader;
class __attribute__((annotate("$clingAutoload$WCLTreeWriter.hh")))  WCLTreeWriter;
class __attribute__((annotate("$clingAutoload$SandBoxPMTcoverage.hh")))  SandBoxPMTcoverage;
class __attribute__((annotate("$clingAutoload$WCSimTrueLight.hh")))  WCSimTrueLight;
class __attribute__((annotate("$clingAutoload$WCSimTruePart.hh")))  WCSimTruePart;
class __attribute__((annotate("$clingAutoload$WCSimTrueCapture.hh")))  WCSimTrueCapture;
class __attribute__((annotate("$clingAutoload$HighEReco.hh")))  HighEReco;
class __attribute__((annotate("$clingAutoload$LowEReco.hh")))  LowEReco;
class __attribute__((annotate("$clingAutoload$LikelihoodGenerator.hh")))  LikelihoodGenerator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "WCLRootDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "WCSimRecoObjectTable.hh"
#include "WCSimRecoDigit.hh"
#include "SandBoxPMTcoverage.hh"
#include "WCSimTrueLight.hh"
#include "WCSimTruePart.hh"
#include "WCSimTrueCapture.hh"
#include "WCLTreeReader.hh"
#include "WCLTreeWriter.hh"
#include "WChRecoLite.hh"
#include "HighEReco.hh"
#include "LowEReco.hh"
#include "LikelihoodGenerator.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"HighEReco", payloadCode, "@",
"LikelihoodGenerator", payloadCode, "@",
"LowEReco", payloadCode, "@",
"SandBoxPMTcoverage", payloadCode, "@",
"WCLTreeReader", payloadCode, "@",
"WCLTreeWriter", payloadCode, "@",
"WCSimRecoDigit", payloadCode, "@",
"WCSimRecoObjectTable", payloadCode, "@",
"WCSimTrueCapture", payloadCode, "@",
"WCSimTrueLight", payloadCode, "@",
"WCSimTruePart", payloadCode, "@",
"WChRecoLite", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("WCLRootDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_WCLRootDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_WCLRootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_WCLRootDict() {
  TriggerDictionaryInitialization_WCLRootDict_Impl();
}
