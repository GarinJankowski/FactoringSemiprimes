from tkinter import *


class ModCalc:
    def __init__(self, modulo):
        self.modulo = modulo

    ENTRY_WIDTH = 7

    def add(self, num1, num2):
        return self.mod(num1 + num2)

    def multiply(self, factor1, factor2):
        return self.mod(factor1 * factor2)

    def divide(self, dividend, divisor):
        return self.mod(dividend / divisor)

    def power(self, num, power):
        return self.mod(num ** power)

    # creates a tkinter table of either addition, multiplication, exponentiation within a mod
    # choose a range for the x and y axis and it operates on them in whatever the ModCalc's mod is
    def create_table(self, symbol, range1, range2, inc1=1, inc2=1):
        root = Tk()
        root.title("mod " + str(self.modulo))

        start1, end1 = range1
        start2, end2 = range2

        entry = Entry(root, width=self.ENTRY_WIDTH)
        entry.grid(row=0, column=0)
        entry.insert(END, "x" + symbol + "y")

        row = 1
        column = 0
        for y in range(start2, end2, inc2):
            entry = Entry(root, width=self.ENTRY_WIDTH)
            entry.grid(row=row, column=0)
            entry.insert(END, y)
            row += 1
        column += 1
        for x in range(start1, end1, inc1):
            entry = Entry(root, width=self.ENTRY_WIDTH)
            entry.grid(row=0, column=column)
            entry.insert(END, x)
            row = 1
            for y in range(start2, end2, inc2):
                if symbol == "+":
                    value = self.add(x, y)
                elif symbol == "*":
                    value = self.multiply(x, y)
                elif symbol == "/":
                    value = self.divide(x, y)
                elif symbol == "**":
                    value = self.power(x, y)
                entry = Entry(root, width=self.ENTRY_WIDTH)
                entry.grid(row=row, column=column)
                entry.insert(END, value)
                row += 1
            column += 1
        root.mainloop()

    # mod 0 means no mod
    def mod(self, number):
        if self.modulo != 0:
            number %= self.modulo
        return number

    def set_mod(self, modulo):
        self.modulo = modulo


def mod(number, modulo):
    if modulo != 0:
        number %= modulo
    return number
