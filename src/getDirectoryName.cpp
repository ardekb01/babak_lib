#include <string.h>

void getDirectoryName(const char *pathname, char *dirname)
{
	int n;
	int i;

	n = (int)strlen(pathname);

	for(i=n-1;i>=0;i--)
	if( pathname[i] == '/') break;

	if(i==-1)
	{
		dirname[0]='.';
		dirname[1]='\0';
	}
	else
	{
		strncpy(dirname, pathname, i);
		dirname[i]='\0';
	}

	return;
}
