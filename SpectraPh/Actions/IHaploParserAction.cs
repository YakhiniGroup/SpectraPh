namespace SpectraPh.Actions
{
    public interface IHaploParserAction
    {
        /// <summary>
        /// Return the usage string which should be displayed on the console
        /// </summary>
        string Usage();

        /// <summary>
        /// The name of the command which is used to invoke this action
        /// </summary>
        string Command();

        /// <summary>
        /// Does all the work, received the arguments from the command line (without the process name & action name)
        /// </summary>
        void Go(string[] args);
    }
}
