CFLAGS = -I../libastro

run-test: test
	./test > output.actual
	diff -u output.expected output.actual

test: test.o ../libastro/libastro.a
	$(CC) $< ../libastro/libastro.a -lm -o $@

test.o: test.c

../libastro/libastro.a: .FORCE
	make -C ../libastro

.PHONY: .FORCE
