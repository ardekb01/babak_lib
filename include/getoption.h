#ifndef GETOPTION_H
#define GETOPTION_H

extern int optind;
extern const char *optArg;

struct CmdOption
{
    const char *name;
    int has_arg;
    int val;
};

int getoption(int argc, char *argv[], const struct CmdOption *options);

#endif
