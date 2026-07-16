#ifndef GETOPTION_H
#define GETOPTION_H

extern int optind;
extern char *optarg;

struct CmdOption
{
    const char *name;
    int has_arg;
    int val;
};

int getoption(int argc, char *argv[], struct CmdOption *options);

#endif
