def display_onscreen_message(message='', color='red'):
    full_message = [message]
    if len(message) > 30:
        full_message = []
        split_message = message.split(' ')
        line_length = 0
        new_line = ''
        e = 0
        for i in range(len(split_message)):
            print(i, len(split_message) - 1)
            line_length += len(split_message[i]) + 1
            if line_length >= 30:
                new_line = ' '.join(split_message[e:i])
                print(new_line, len(new_line))
                full_message.append(new_line)
                new_line = ''
                line_length = 0
                e = i
                print(split_message[e:])
            if int(i) == int(len(split_message) - 1):
                print('k')
                new_line = ' '.join(split_message[e:])
                print(new_line, len(new_line))
                full_message.append(new_line)
            print int(i), int(len(split_message) - 1)
    for line in full_message:
        print(line)


if __name__ == '__main__':
    display_onscreen_message(message='There is still no wavelength solution!')
