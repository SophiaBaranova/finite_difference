void error(int code)
{
	switch (code)
	{
	case 1:
		printf("file opening error\n");
		exit(1);
	case 2:
		printf("memory allocation error\n");
		exit(2);
	case 3:
		printf("invalid calculations\n");
		exit(4);
	default:
		exit(code);
	}
}
