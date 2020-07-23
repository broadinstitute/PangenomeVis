/**
 * Strings Splitter wrapper.
 *
 * @author Htmlstream
 * @version 1.0
 *
 */
;(function ($) {
	'use strict';
	$.HSCore.components.HSStringsUnifier = {
		/**
		 *
		 *
		 * @var Object _baseConfig
		 */
		_baseConfig: {
			splitter: ", "
		},

		/**
		 *
		 *
		 * @var jQuery pageCollection
		 */
		pageCollection: $(),

		/**
		 * Initialization of Strings Splitter wrapper.
		 *
		 * @param String selector (optional)
		 * @param Object config (optional)
		 *
		 * @return jQuery pageCollection - collection of initialized items.
		 */

		init: function (selector, config) {
			this.collection = selector && $(selector).length ? $(selector) : $();
			if (!$(selector).length) return;

			this.config = config && $.isPlainObject(config) ?
				$.extend({}, this._baseConfig, config) : this._baseConfig;

			this.config.itemSelector = selector;

			this.initStringsUnifier();

			return this.pageCollection;
		},

		initStringsUnifier: function () {
			//Variables
			var $self = this,
				config = $self.config,
				collection = $self.pageCollection;

			//Actions
			this.collection.each(function (i, el) {
				//Variables
				var $this = $(el),
					isInput = $this[0].attributes.value,

					options = JSON.parse(el.getAttribute('data-su-options')),

					resultString = '';

				$(options.data).each(function (i, el) {

					var targetValues = 0;

					$(el.target).each(function () {

						targetValues += parseInt($(this)[0].attributes.value ? $(this).val() : $(this).text());

					});

					if (targetValues !== 0) {

						resultString += i === 0 ? '' : (options.splitter ? options.splitter : config.splitter);
						resultString += targetValues;
						resultString += ' ';
						resultString += el.title ? el.title : '';
						resultString += el.titleMultiplier && targetValues !== 1 ? el.titleMultiplier : '';

					}

				});

				if (isInput) {

					$this.val(resultString);

				} else {

					$this.text(resultString);

				}

				//Actions
				collection = collection.add($this);
			});
		}
	};
})(jQuery);
