#include <iostream>
#include <pthread.h>

void *hello(void *args) {
	char *mesg = (char *) args;
	const int id = pthread_self();
	fprintf(stdout, "[thread=%02d][ LOG ] Hello with message: %s\n", id, mesg);
	return NULL;
}

int main() {

	pthread_t thread1, thread2;

	pthread_create(&thread1, NULL, &hello, (void *) "from one");
	pthread_create(&thread2, NULL, &hello, (void *) "from two");

	pthread_join(thread1, NULL);
	pthread_join(thread2, NULL);

	return EXIT_SUCCESS;
}
