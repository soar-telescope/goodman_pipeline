import PySimpleGUI as sg


def display_help_text():
    sg.theme('DarkGrey13')
    layout = [[[sg.Text("Help:")],
               [sg.MLine(size=(50, 10), key='-HELP-')],
               [sg.Button("Ok")]]]

    window = sg.Window("Help", layout)
    while True:
        event, values = window.read()
        window['-HELP-'].print("Ctrl+q: Exit the interactive mode.")

        if event == sg.WIN_CLOSED or event == 'Ok': # if user closes window or clicks cancel
            break


    window.close()


if __name__ == '__main__':
    display_help_text()
