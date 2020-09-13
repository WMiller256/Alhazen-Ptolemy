#Name of program
MAIN     = alhazen-ptolemy

ABS      = ./
BIN      = ~/bin
BUILD    = ./
RM       = /bin/rm -f
MV       = /bin/mv -f
CFLAGS   = -std=c++17 -O3
CC       = /usr/bin/c++ $(CFLAGS)

LIBS     =
                    
LFLAGS   = -Wall -Wl,-rpath,/usr/local/lib
LIBDIRS  = $(LFLAGS) -L/usr/local/lib/

#Output coloring
GREEN   = \033[1;32m
CYAN    = \033[36m
BLUE    = \033[1;34m
BRIGHT  = \033[1;37m
WHITE   = \033[0;m
MAGENTA = \033[35m
YELLOW  = \033[33m
RED     = \033[91m

#Source files
OBJS   = $(BUILD)/Alhazen-Ptolemy.o

#Builds
$(BUILD)/%.o: %.c++
	@printf "[$(CYAN)Building$(WHITE)]   $(BRIGHT)$<$(WHITE) - $(MAGENTA)Object$(WHITE)\n"
	cd $(ABS); $(CC) -c -o $@ $<
	@printf "[$(GREEN) Built  $(WHITE)]   $(BRIGHT)$<$(WHITE) - $(MAGENTA)Object$(WHITE)\n"

Alhazen-Ptolemy.o: Alhazen-Ptolemy.c++
	@printf "[$(CYAN)Building$(WHITE)]   $(BRIGHT)$<$(WHITE) - $(MAGENTA)Object$(WHITE)\n"
	cd $(ABS); $(CC) Alhazen-Ptolemy.c++ -o $(BIN)/$(MAIN)
	@printf "[$(GREEN) Built  $(WHITE)]   $(BRIGHT)$<$(WHITE) - $(MAGENTA)Object$(WHITE)\n"

clean:
	$(RM) *.core $(BUILD)/*.o *.d *.stackdump

#Disable command echoing, reenabled with make verbose=1
ifndef verbose
.SILENT:
endif
