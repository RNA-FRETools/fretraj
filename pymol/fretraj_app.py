#!/usr/bin/python3

import sys

if sys.version_info.major < 3:
    import Tkinter as tk
else:
    import tkinter as tk


class App(tk.Frame):
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.pack()
        self.master.title("Hello World")
        self.master.resizable(False, False)
        self.master.tk_setPalette(background='#ececec')

        self.master.protocol('WM_DELETE_WINDOW', self.click_cancel)
        self.master.bind('<Return>', self.click_ok)
        self.master.bind('<Escape>', self.click_cancel)

        x = int((self.master.winfo_screenwidth() - self.master.winfo_reqwidth()) / 2)
        y = int((self.master.winfo_screenheight() - self.master.winfo_reqheight()) / 3)
        self.master.geometry("+{:d}+{:d}".format(x, y))

        self.master.config(menu=tk.Menu(self.master))

        dialog_frame = tk.Frame(self)
        dialog_frame.pack(padx=20, pady=15)

        tk.Label(dialog_frame, text="This is your first GUI. (highfive)").pack()

        button_frame = tk.Frame(self)
        button_frame.pack(padx=15, pady=(0, 15), anchor='e')

        tk.Button(button_frame, text='OK', default='active', command=self.click_ok).pack(side='right')

        tk.Button(button_frame, text='Cancel', command=self.click_cancel).pack(side='right')

    def click_ok(self, event=None):
        print("The user clicked 'OK'")

    def click_cancel(self, event=None):
        print("The user clicked 'Cancel'")
        self.master.destroy()


if __name__ == '__main__':
    root = tk.Tk()
    root.lift()
    app = App(root)
    app.mainloop()
