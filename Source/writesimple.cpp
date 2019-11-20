/* Write output in simple format. */

#include "write.h"
#include <string.h> /* FACTS */
#include <time.h>

#define F_FLOAT "%.25f"

/* --------------- Simple ----------------- */
void Simple(void) {
  FILE *OutFile;

  if (ExistParameter('=')) {
    if (NullName(NameParameter('=')))
      return;
    OutFile = OpenOutput(NameParameter('='));
  }
  // fprintf(OutFile, "BASE " F_FLOAT "\n", Sn);

  /* --------------------- AC bus results -----------------------------*/
  // fprintf(OutFile, "ACBUS %d\n", Nac);
  for (auto ACptr = dataPtr->ACbus; ACptr != nullptr; ACptr = ACptr->Next) {
    fprintf(OutFile, "%6d %12s ", ACptr->Num, ACptr->Name);

    auto delta = ACptr->Ang;
    auto vals = delta >= 0 ? 1.0 : -1.0;
    if (fabs(delta) > 2 * PI)
      delta = delta - vals * floor(fabs(delta) / (2 * PI)) * 2 * PI;
    if (fabs(delta) > PI)
      delta = delta - vals * 2 * PI;
    ACptr->Ang = delta;
    fprintf(OutFile, F_FLOAT " " F_FLOAT " ", ACptr->V, ACptr->Ang * (180.0 / PI));

    auto Pl = (ACptr->Pn + lambda * ACptr->Pnl) * pow(ACptr->V, ACptr->a) + (ACptr->Pz + lambda * ACptr->Pzl) * ACptr->V * ACptr->V;
    auto Ql = (ACptr->Qn + lambda * ACptr->Qnl) * pow(ACptr->V, ACptr->b) + (ACptr->Qz + lambda * ACptr->Qzl) * ACptr->V * ACptr->V;
    fprintf(OutFile, F_FLOAT " " F_FLOAT " ", Pl * Sn, Ql * Sn);

    fprintf(OutFile, F_FLOAT " " F_FLOAT "\n", ACptr->PG * Sn, ACptr->Qg * Sn);
  }
  fclose(OutFile);
}
