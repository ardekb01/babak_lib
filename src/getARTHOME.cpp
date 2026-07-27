#include <cstdlib>

#define GETARTHOME

const char *ARTHOME = nullptr;

// Get the value of the ARTHOME environment variable.
// The getenv() function searches the environment list for a string
// that matches "ARTHOME". It returns a pointer to the value in the
// environment, or nullptr if there is no match.
bool getARTHOME()
{
   ARTHOME = getenv("ARTHOME");

   if(ARTHOME == nullptr)
   {
      return false;
   }

   return true;
}
